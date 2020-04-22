%voltage sources
gnd = vsrc('gnd',0,1);
vdd = vsrc('vdd',1,1);
clk = vtanhClock('clk',1e10,-pi/2,0.5,0.5,4,3e-10,0,1);
clkbar = vtanhClock('clkbar',1e10,pi/2,0.5,0.5,4,3e-10,1,1);
din = vtanhDelay('din',0.5,0.5,5e10,true);

%simulation parameters
timeSpan = [0 1.35e-9];
tCrit = 1.3e-9;

sync = TWO_FLOP_PG_SYNC('pgSync_test',1,1);

mask = [synch.cctPath.PGFF_0.PGL_0.q,...
        synch.cctPath.PGFF_0.q,...
        synch.cctPath.PGFF_1.PGL_0.q,...
        synch.cctPath.q];

numTx = 42;
numTxRobust = 78;

nomWidth = 450e-7;
minWidth = 200e-7;

lastLatch = 33:1:numTx;

uVcrit = zeros(1,synch.nodeNum);
uNodes = [synch.cctPath.PGFF_1.PGL_0.q,...
          synch.cctPath.PGFF_1.PGL_0.x0];
uVcrit(uNodes) = [1,-1];

transistorNumbering = 1:1:numTx;
nDevices = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41];
pDevices = [2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42];
transistorWidths = ones(numTx,1)*nomWidth;
scale = (sum(ones(1,numTxRobust)*minWidth*4))/(sum(transistorWidths));
transistorWidths = scale*transistorWidths;
wMin = ones(numTx,1)*minWidth;
wMin(lastLatch) = nomWidth*0.999;
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
while(norm(transistorWidths-transistorWidthPrev) > norm(wStart)/numTx*tol)
    close all
    sync = TWO_FLOP_PG_SYNC('pgSync_test',transistorWidths,1);
    
    nb = nestedBisection(sync,{sync.vdd,sync.gnd,sync.clk,sync.clkbar,sync.d},...
        {vdd,gnd,clk,clkbar,vdd},[5.5e-10 3e-10],{sync.d,din},clk,mask,'EKV','PTM 45nmHP');
    
    dataSync = nb.simulate(strcat('twoFlop_PG_Sync_iter',num2str(i)),timeSpan,tCrit);
    
    nbAnalysis = nestedBisectionAnalysis(sync,{sync.vdd,sync.gnd,sync.clk,sync.clkbar,sync.d},...
        {vdd,gnd,clk,clkbar,vdd},{sync.d,din},clk,mask(end-1),'EKV','PTM 45nmHP');
    
    nbAnalysis.bisectionAnalysis(strcat('twoFlop_PG_SyncAnalysis_',num2str(i)),dataSync,uVcrit,tCrit);
    
    load(strcat('twoFlop_PG_SyncAnalysis_',num2str(i),'.mat'))
    teola = t(end);
    saveTeola(i) = teola;
    fprintf('Current tCrit: %d\n',saveTeola(i));
    dGdwRaw = nbAnalysis.dGdw(splMeta,Jac_t,dhda_t,dataTransition,beta_t,[t(1),teola],uVcrit);
    
    dGdw = zeros(numTx,1);
    dGdw(nDevices) = dGdwRaw(1:length(nDevices));
    dGdw(pDevices) = dGdwRaw(length(nDevices)+1:end);
    
    dGdw = -1/g*dGdw;
    %
    dGdw(lastLatch) = 0;
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
    
    i = i +1;
    txWidthsSave(:,i) = transistorWidths;
    gradientSave(:,i) = dGds;
    fprintf('Error ||txCurrent - txPrev|| %d\n',norm(transistorWidths-transistorWidthPrev))
end

ind = nnz(saveTeola);
tCrits = saveTeola(1:ind);
txWidths = txWidthsSave(:,1:ind);
save('twoFlopSynchronizerTxWidths','saveTeola','txWidths','gradientSave','gradientRawSave');

    
