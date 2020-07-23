
%voltage sources
gnd = vsrc('gnd',0,1);
vdd = vsrc('vdd',1,1);
clkbar = vtanhClock('clk',1e10,-pi/2,0.5,0.5,4,3e-10,0,1);
clk = vtanhClock('clkbar',1e10,pi/2,0.5,0.5,4,3e-10,1,1);
din = vtanhDelay('din',-0.5,0.5,5e10,true);
reset = vpulse('reset',[0,1],[0.1e-12,1e-10,0.1e-12,1.5e-9],1e-10);
%simulation parameters
timeSpan = [0 1.5e-9];
tCrit = 1.3e-9;

sync = TWO_FLOP_ROBUST_JAMB_SYNC('robustJambSync_test',1,1);
mask = [sync.cctPath.RJFF_0.RJL_0.y,sync.cctPath.RJFF_0.RJL_1.y,...
        sync.cctPath.RJFF_1.RJL_0.y,sync.cctPath.RJFF_1.RJL_1.y];

nomWidth = 450e-7;
minWidth = 200e-7;
numTx = 78;
margin = 0.2;

%%% Transistor sizing starting point (near optimal point)
%Buffer
invBuf = [1150e-7;250e-7];
%Latch
invCCfwd  = [1400e-7;300e-7];
invCCBck  = [1400e-7;250e-7];
invOut    = [275e-7;1140e-7];
nand0 = [250e-7;1140e-7;1140e-7;280e-7];
nDin        = 330e-7;
nClk        = 340e-7;
nReset      = 250e-7;
nFilter0    = 250e-7;
nFilter1    = 250e-7;
pBoost0     = 215e-7;
pBoost1     = 215e-7;
pFilter0    = 260e-7;
pFilter1    = 260e-7;

latchSize = [invCCfwd;invCCBck;invOut;nand0;nDin;nClk;nReset;nFilter0;nFilter1;pBoost0;pBoost1;pFilter0;pFilter1];
ff = [latchSize;latchSize];
transistorWidths = [invBuf;ff;ff];


%%%%

uVcrit = zeros(1,39);
uVcrit(sync.cctPath.RJFF_1.RJL_0.x) = 1;
uVcrit(sync.cctPath.RJFF_1.RJL_0.y) = -1;


nDevices   = [1, 3,5,7,9,10,13,14,15,16,17, 22,24,26,28,29,32,33,34,35,36,  41,43,45,47,48,51,52,53,54,55,  60,62,64,66,67,70,71,72,73,74]; 
pDevices  =  [2, 4,6,8,11,12,18,19,20,21,   23,25,27,30,31,37,38,39,40,     42,44,46,49,50,56,57,58,59,     61,63,65,68,69,75,76,77,78];

crossCouplePairs = [3,4,5,6,22,23,24,25,41,42,43,44,60,61,62,63];
pulldowns        = [13,14,15, 32,33,34, 51,52,53, 67,70,71];
resetTxs         = [15,34,53,71];

lastLatch = 60:numTx;

v0 = zeros(1,39);
v0(mask) = 1.0;

%transistorWidths(pulldowns) = minWidth;
scale = (sum(ones(1,numTx)*minWidth*4))/sum(transistorWidths);
transistorWidths = scale*transistorWidths;

transistorWidths = pullDownSizing(crossCouplePairs,pulldowns,transistorWidths,margin);

pulldowns = pulldowns(1:9);
crossCouplePairs = crossCouplePairs(1:12);

wMin = ones(numTx,1)*minWidth;
wMin(lastLatch) = transistorWidths(lastLatch)*0.999;

f_s = @(txWidth) txWidth./wMin - 1./(txWidth./wMin);
finv_s = @(s) wMin.*(1+(s+sqrt(s.^2+4))/2);

u = transistorWidths./wMin - 1;
s = u - 1./u;
sPrev = zeros(numTx,1);

wStart = transistorWidths;
totalWidth = sum(ones(1,numTx)*minWidth*4);
transistorWidthPrev = zeros(numTx,1);
dGdsPrev = zeros(numTx,1);
tol = 0.5;
eps = 0.01*minWidth;
i = 1;
saveTeola = zeros(1,1000);

txWidthsSave = zeros(length(transistorWidths),1000);
txWidthsSave(:,i) = transistorWidths;
gradientSave = zeros(length(transistorWidths),1000);
gradientRawSave = zeros(length(transistorWidths),1000);

while(norm(transistorWidths-transistorWidthPrev) > norm(wStart)/numTx*tol || i < 3)
    close all
    sync = TWO_FLOP_ROBUST_JAMB_SYNC('robustJambSync_test',transistorWidths,1);
    
    nb = nestedBisection(sync,{sync.vdd,sync.gnd,sync.clk,sync.clkbar,sync.d,sync.reset},...
        {vdd,gnd,clk,clkbar,vdd,reset},[6e-10 2e-10],{sync.d,din},clk,mask,'EKV','PTM 45nmHP');
    
    dataOneLatch = nb.simulate(strcat('twoFlop_RJL_Sync_iter',num2str(i)),timeSpan,tCrit,v0');
    
    nbAnalysis = nestedBisectionAnalysis(sync,{sync.vdd,sync.gnd,sync.clk,sync.clkbar,sync.d,sync.reset},...
        {vdd,gnd,clk,clkbar,vdd,reset},{sync.d,din},clk,mask(end),'EKV','PTM 45nmHP');
    
    nbAnalysis.bisectionAnalysis(strcat('twoFlop_RJL_SyncAnalysis_',num2str(i)),dataOneLatch,uVcrit,tCrit);
    
    load(strcat('twoFlop_RJL_SyncAnalysis_',num2str(i),'.mat'))
    teola = t(end);
    saveTeola(i) = teola;
    fprintf('Current tCrit: %d\n',saveTeola(i));
    dGdwRaw = nbAnalysis.dGdw(splMeta,Jac_t,dhda_t,dataTransition,u_t,beta_t,g_t,[t(1),teola],t,strcat('gradientData_',num2str(i)));
    
    dGdw = zeros(numTx,1);
    dGdw(nDevices) = dGdwRaw(1:length(nDevices));
    dGdw(pDevices) = dGdwRaw(length(nDevices)+1:end);
    
    logdGdw = -1/g_t(teola)*dGdw;
    
    logdGdw(lastLatch) = 0;
    gradientRawSave(:,i) = logdGdw;
    
    r = ones(size(transistorWidths))/sqrt(length(transistorWidths));
    dgdwFix = logdGdw - r.*(r.*logdGdw);
    dwds = wMin.*(1+s./sqrt(s.^2+4))/2;
    dGds = dgdwFix.*dwds;
    
    gamma = abs((s - sPrev)'*(dGds - dGdsPrev))/(norm(dGds - dGdsPrev)^2);
    
    transistorWidthPrev = transistorWidths;
    dGdsPrev = dGds;
    sPrev = s;
    
    s = s - 0.1*gamma*dGds;
    s(lastLatch) = -1e10;
    widBudget = @(sfix) sum(finv_s(s + sfix*ones(size(s)))) - totalWidth;
    sfix = fzero(widBudget,1);
    
    s = s + sfix*ones(size(s));
    
    transistorWidths = finv_s(s);
    
    transistorWidths(1:59) = pullDownSizing(crossCouplePairs,pulldowns,transistorWidths(1:59),margin);
    
    s = f_s(transistorWidths);
    
    save(strcat('gradientState_',num2str(i)),'transistorWidths','transistorWidthPrev','s','sPrev','dGds','dGdwRaw','dGdw','gamma');

    gradientSave(:,i) = dGds;
    i = i +1;
    txWidthsSave(:,i) = transistorWidths;
    fprintf('Error ||txCurrent - txPrev|| %d\n',norm(transistorWidths-transistorWidthPrev))
end

ind = nnz(saveTeola);
tCrits = saveTeola(1:ind);
txWidths = txWidthsSave(:,1:ind);

save('twoFlopSynchronizerTxWidths','tCrits','txWidths');
