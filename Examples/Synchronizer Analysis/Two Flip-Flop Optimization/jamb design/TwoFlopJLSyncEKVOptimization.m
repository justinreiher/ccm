% Configure the testbench
margin = 0.2;

%voltage sources
gnd = vsrc('gnd',0,1);
vdd = vsrc('vdd',1,1);
clkbar = vtanhClock('clk',1e10,-pi/2,0.5,0.5,4,3e-10,0,1);
clk = vtanhClock('clkbar',1e10,pi/2,0.5,0.5,4,3e-10,1,1);
din = vtanhDelay('din',-0.5,0.5,5e10,true);
reset = vpulse('reset',[0,1],[0.1e-12,1e-10,0.1e-12,1.5e-9],1e-10);

sync = TWO_FLOP_JAMB_SYNC('jambSync_test',1,1);

%simulation parameters
timeSpan = [0 1.35e-9];
tCrit = 1.3e-9;
mask = [sync.cctPath.JFF_0.JL_0.q,sync.cctPath.JFF_0.q,...
        sync.cctPath.JFF_1.JL_0.q,sync.cctPath.q];

nomWidth = 450e-7;
pullDownFactor = 1.2;
numTx = 32;
numTxRobust = 78;
minWidth = 200e-7;

lastLatch = 30:1:numTx;
pullDownsEnd = [30,31,32];

v0 = zeros(1,20);
v0(mask) = 1.0;

uVcrit = zeros(1,20);
uVcrit(sync.cctPath.JFF_1.JL_0.q) = -1;
uVcrit(sync.cctPath.JFF_1.JL_0.qbar) = 1;
transistorNumbering = 1:1:numTx;
nDevices = [1,3, 5,7,9,10,11, 12,14,16,17,18, 19,21,23,24,25, 26,28,30,31,32];
pDevices = [2,4, 6,8,         13,15,          20,22,          27,29];
transistorWidths = ones(numTx,1)*nomWidth;
scale = (sum(ones(1,numTxRobust)*minWidth*4))/(sum(transistorWidths));
transistorWidths = scale*transistorWidths;
pulldowns = [9,10,11,16,17,18,23,24,25,30,31,32];
resetTxs  = [11,18,25,32];
lastPullDowns = [30,31,32];
crossCouplePairs = [5,6,7,8,12,13,14,15,19,20,21,22,26,27,28,29];

transistorWidths = pullDownSizing(crossCouplePairs,pulldowns,transistorWidths,margin);

wMin = ones(numTx,1)*minWidth;
wMin(lastLatch) = nomWidth*0.999;
wMin(lastPullDowns) = nomWidth*pullDownFactor*0.999;

f_s = @(txWidth) txWidth./wMin - 1./(txWidth./wMin);
finv_s = @(s) wMin.*(1+(s+sqrt(s.^2+4))/2);

u = transistorWidths./wMin - 1;
s = u - 1./u;
sPrev = zeros(numTx,1);

wStart = transistorWidths;
totalWidth = sum(transistorWidths);
transistorWidthPrev = zeros(numTx,1);
dGdsPrev = zeros(numTx,1);
tol = 0.2;
eps = 0.01*minWidth;
i = 1;
saveTeola = zeros(1,1000);

txWidthsSave = zeros(length(transistorWidths),1000);
txWidthsSave(:,i) = transistorWidths;
gradientSave = zeros(length(transistorWidths),1000);
gradientRawSave = zeros(length(transistorWidths),1000);

while(norm(transistorWidths-transistorWidthPrev) >= norm(wStart)*tol/numTx)
    close all

    sync = TWO_FLOP_JAMB_SYNC('jambSync_test',transistorWidths,1);
    
     nb = nestedBisection(sync,{sync.vdd,sync.gnd,sync.clk,sync.clkbar,sync.d,sync.reset},...
         {vdd,gnd,clk,clkbar,vdd,reset},[5.5e-10 2e-10],{sync.d,din},clk,mask,'EKV','PTM 45nmHP');
     
     dataSynch = nb.simulate(strcat('twoFlop_jambSync_iter',num2str(i)),timeSpan,tCrit,v0');
     
     nbAnalysis = nestedBisectionAnalysis(sync,{sync.vdd,sync.gnd,sync.clk,sync.clkbar,sync.d,sync.reset},...
         {vdd,gnd,clk,clkbar,vdd,reset},{sync.d,din},clk,mask(end-1),'EKV','PTM 45nmHP');
     
     nbAnalysis.bisectionAnalysis(strcat('twoFlop_jambSyncAnalysis_',num2str(i)),dataSynch,uVcrit,tCrit);
    
    load(strcat('twoFlop_jambSyncAnalysis_',num2str(i),'.mat'))
    teola = t(end);
    gSign = sign(g);
    saveTeola(i) = teola;
    fprintf('Current tCrit: %d\n',saveTeola(i));
    dGdwRaw = nbAnalysis.dGdw(splMeta,Jac_t,dhda_t,dataTransition,beta_t,[t(1),teola],uVcrit);
    
    dGdw = zeros(numTx,1);
    dGdw(nDevices) = dGdwRaw(1:length(nDevices));
    dGdw(pDevices) = dGdwRaw(length(nDevices)+1:end);
    
    dGdw = -1/g*dGdw;
    %
    dGdw(lastLatch) = 0;
    dGdw(resetTxs) = 0;
    gradientRawSave(:,i) = dGdw;
    
    r = ones(size(transistorWidths))/sqrt(length(transistorWidths));
    dgdwFix = dGdw - r.*(r.*dGdw);
    dwds = wMin.*(1+s./sqrt(s.^2+4))/2;
    dGds = dgdwFix.*dwds;
    
    gamma = abs((s - sPrev)'*(dGds - dGdsPrev))/(norm(dGds - dGdsPrev)^2);
    
    transistorWidthPrev = transistorWidths;
    dGdsPrev = dGds;
    sPrev = s;
    
    s = s - gamma*dGds;
    
    widBudget = @(sfix) sum(finv_s(s + sfix*ones(size(s)))) - totalWidth;
    sfix = fzero(widBudget,1);
    
    transistorWidths = finv_s(s+sfix*ones(size(s)));
    
     transistorWidths = pullDownSizing(crossCouplePairs,pulldowns,transistorWidths,margin);
    
    save(strcat('gradientState_',num2str(i)),'transistorWidths','transistorWidthPrev','s','sPrev','dGds','dGdwRaw','dGdw','gamma');

    gradientSave(:,i) = dGds;
    i = i +1;
    txWidthsSave(:,i) = transistorWidths;
    fprintf('Error ||txCurrent - txPrev|| %d\n',norm(transistorWidths-transistorWidthPrev))
end

ind = nnz(saveTeola);
tEolas = saveTeola(1:ind);
txWidths = txWidthsSave(:,1:ind);

save('twoFlopSynchronizerTxWidths','saveTeola','txWidths','gradientSave','gradientRawSave');




