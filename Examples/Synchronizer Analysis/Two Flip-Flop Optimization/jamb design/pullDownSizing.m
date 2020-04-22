function transistorWidths = pullDownSizing(crossCoupledWidths,pullDowns,transistorWidths,margin)

txNumbering = 1:length(transistorWidths);

totalWidths = sum(transistorWidths);
ccTxSizes   = transistorWidths(crossCoupledWidths);
pullDownTxSizes = transistorWidths(pullDowns);
excludePullDownInd = setdiff(txNumbering,pullDowns);

timeSpan = [0,0.1e-9];

slope = 1/timeSpan(2);

%voltage sources
gnd = vsrc('gnd',0,1);
vdd = vsrc('vdd',1,1);
din = vsrcLinear('din',slope,true);

margin = 1+margin;


pullDownWidths = zeros(size(pullDowns));


for i = 1:length(crossCoupledWidths)/4
    
    wid = ccTxSizes(1+(i-1)*4:i*4);
    widPullDowns = pullDownTxSizes(1 +(i-1)*3:i*3);
    wid0 = wid(1:2); %forward CC inverter
    wid1 = wid(3:4); %backward CC inverter
          
    n_params0 = ekv_model_params_45nmHP('n');
    n_params1 = ekv_model_params_45nmHP('n');
    n_params0.W = wid0(1);
    n_params1.W = wid1(1);
    
    p_params0 = ekv_model_params_45nmHP('p');
    p_params1 = ekv_model_params_45nmHP('p');
    p_params0.W = wid0(2);
    p_params1.W = wid1(2);
    
    inv0 = INV4('tf',wid0,1);
    inv1 = INV4('tf',wid1,1);
    
    tfInv0Bench = testbench(inv0,{inv0.i,inv0.vdd,inv0.gnd},{din,vdd,gnd},'EKV','PTM 45nmHP');
    tfInv1Bench = testbench(inv1,{inv1.i,inv1.vdd,inv1.gnd},{din,vdd,gnd},'EKV','PTM 45nmHP');
    
    [ti0,vi0] = tfInv0Bench.simulate(timeSpan,[0;0;0;0.99]);
    [ti1,vi1] = tfInv1Bench.simulate(timeSpan,[0;0;0;0.99]);
    
    tf0 = spline(din.V(ti0),vi0(:,4));
    tf1 = spline(din.V(ti1),vi1(:,4));
    
    tf0_v = @(v) ppval(tf0,v);
    tf1_v = @(v) ppval(tf1,v);
    
    Vx = linspace(0,1);
    Vy = linspace(0,1);
    
    ix = -ekv_Ids(p_params1,[Vx;tf0_v(Vx);ones(1,100);ones(1,100)]) ...
          - ekv_Ids(n_params1,[Vx;tf0_v(Vx);zeros(1,100);zeros(1,100)]);
      
    iy = - ekv_Ids(p_params0,[Vy;tf1_v(Vx);ones(1,100);ones(1,100)]) ...
          - ekv_Ids(n_params0,[Vy;tf1_v(Vx);zeros(1,100);zeros(1,100)]);
    
    [ix,ind] = max(ix);
    v0Max = Vx(ind);
    
    
    nmosNom = ekv_model_params_45nmHP('n');
    nmosD = nmosNom;
    nmosClk = nmosNom;
    
    nmosD.W = widPullDowns(1);
    nmosClk.W = widPullDowns(2);
    
    iz = @(z) ekv_Ids(nmosD,[v0Max;1.0;z;0]) - ekv_Ids(nmosClk,[z;1.0;0;0]);
    
    vZ = fzero(iz,v0Max);
    
    iTx = ekv_Ids(nmosClk,[vZ;1.0;0;0]);
    
    wScale = ix*margin/iTx;
    
    if(wScale > 1)
        wClk = nmosClk.W*wScale;
        wD   = nmosD.W*wScale;
    else 
        wClk = nmosClk.W;
        wD   = nmosD.W;
    end
    
    nmosCheck = nmosNom;
    nmosCheck.W = wD;
   
    iCheck = ekv_Ids(nmosCheck,[v0Max;1.0;vZ;0]);
    assert(iCheck > ix);
    nmosCheck.W = wClk;
    iCheck = ekv_Ids(nmosCheck,[vZ;1.0;0;0]);
    assert(iCheck > ix); 
    
    pullDownWidths(1+(i-1)*3) = wD;
    pullDownWidths(2+(i-1)*3) = wClk; 
    
    [iy,ind] = max(iy);
    v1Max = Vy(ind);
    nmosR = ekv_model_params_45nmHP('n');
    nmosR.W = widPullDowns(3);
    
    iTx = ekv_Ids(nmosR,[v1Max;1.0;0;0]);
    wScale = iy*margin/iTx;
    
    if(wScale > 1)
        wR = nmosR.W*wScale;
    else
        wR = nmosR.W;
    end
    
    nmosR.W = wR;

    iCheck = ekv_Ids(nmosR,[v1Max;1.0;0;0]);
    assert(iCheck > iy);
    
    pullDownWidths(3+(i-1)*3) = wR;
end

aPullDowns = pullDownWidths'./pullDownTxSizes;

scaleTxs = (totalWidths - sum(pullDownTxSizes.*aPullDowns))/sum(transistorWidths(excludePullDownInd));
transistorWidths(pullDowns) = aPullDowns.*pullDownTxSizes;
transistorWidths(excludePullDownInd) = scaleTxs'.*transistorWidths(excludePullDownInd);
    
