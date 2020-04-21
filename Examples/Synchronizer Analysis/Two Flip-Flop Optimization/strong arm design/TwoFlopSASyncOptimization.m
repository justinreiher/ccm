%voltage sources
gnd = vsrc('gnd',0,1);
vdd = vsrc('vdd',1,1);
clk = vtanhClock('clk',1e10,-pi/2,0.5,0.5,4,3e-10,0,1);
din = vtanhDelay('din',0.5,0.5,5e10,true);

%simulation parameters
timeSpan = [0 1.5e-9];
tCrit = 1.3e-9;

sync = TWO_FLOP_STRONG_ARM_SYNC('saSync_test',1,1);
mask = [sync.cctPath.saFF_0.s,sync.cctPath.saFF_0.r,...
        sync.cctPath.saFF_0.q,sync.cctPath.saFF_1.r,...
        sync.cctPath.q];

%mask = [6,16,15,13,17];
nomWidth = 450e-7;
minWidth = 200e-7; 
numTx = 40;

uVcrit = zeros(1,24);
uVcrit(13) = 1;
uVcrit(14) = -1;
transistorNumbering = 1:1:numTx;
nDevices = [1,3, 5,7, 9,10,11,12, 15,16,19,20   23,25, 27,28,29,30, 33,34,37,38];
pDevices = [2,4, 6,8, 13,14,      17,18,21,22,  24,26, 31,32,       35,36,39,40];

clkPullDowns = [11,29];
preChargeTx = [13,14,31,32];
nandGndTx   = [16,33];
bridgeTx    = [12,30];
transistorWidths = ones(numTx,1)*nomWidth;
%transistorWidths(clkPullDowns) = nomWidth*1.5;
%transistorWidths(nandGndTx)    = nomWidth*1.5;
%transistorWidths(preChargeTx) = minWidth*1.1;
%transistorWidths(bridgeTx) = nomWidth/2;
scale = sum(ones(1,numTx)*minWidth*4)/(sum(transistorWidths));
transistorWidths = transistorWidths*scale;

wMin = ones(numTx,1)*minWidth;



f_s = @(txWidth) txWidth./wMin - 1./(txWidth./wMin);
finv_s = @(s) wMin.*(1+(s+sqrt(s.^2+4))/2);

u = transistorWidths./wMin - 1;
s = u - 1./u;
sPrev = zeros(numTx,1);

wStart = transistorWidths;
totalWidth = sum(transistorWidths);
transistorWidthPrev = zeros(numTx,1);
dGdsPrev = zeros(numTx,1);
tol = 0.1;
eps = 0.01*minWidth;
i = 1;
saveTeola = zeros(1,1000);

txWidthsSave = zeros(length(transistorWidths),1000);
txWidthsSave(:,i) = transistorWidths;
while(norm(transistorWidths-transistorWidthPrev) >= norm(wStart)/(numTx)*tol)
    close all

    sync = TWO_FLOP_STRONG_ARM_SYNC('sARMS ',transistorWidths,1);
    
    nb = nestedBisection(sync,{sync.vdd,sync.gnd,sync.clk,sync.d},...
        {vdd,gnd,clk,vdd},[5.5e-10 2e-10],{sync.d,din},clk,mask,'EKV','PTM 45nmHP');
    
    dataOneLatch = nb.simulate(strcat('twoFlop_SA_Sync_iter',num2str(i)),timeSpan,tCrit);
    
    nbAnalysis = nestedBisectionAnalysis(sync,{sync.vdd,sync.gnd,sync.clk,sync.d},...
        {vdd,gnd,clk,vdd},{sync.d,din},clk,mask(end-1),'EKV','PTM 45nmHP');
    
    nbAnalysis.bisectionAnalysis(strcat('twoFlop_SA_SyncAnalysis_',num2str(i)),dataOneLatch,uVcrit,tCrit);
    
    load(strcat('twoFlop_SA_SyncAnalysis_',num2str(i),'.mat'))
    teola = t(end);
    saveTeola(i) = teola;
    gSign = sign(g);
    if(gSign > 0 )
        gSign = -1;
    else
        gSign = 1;
    end
    fprintf('Current tCrit: %d\n',saveTeola(i));
    dGdwRaw = nbAnalysis.dGdw(splMeta,Jac_t,dhda_t,dataTransition,beta_t,[t(1),tCrit],uVcrit);
    
    dGdw = zeros(numTx,1);
    dGdw(nDevices) = dGdwRaw(1:length(nDevices));
    dGdw(pDevices) = dGdwRaw(length(nDevices)+1:end);
    
%     signGrad = sign(dGdw);
%     log_dGdw = log(abs(dGdw));
%     dGdw = -gSign*signGrad.*log_dGdw;
    dGdw = 1/g*dGdw;
    
    r = ones(size(transistorWidths))/sqrt(length(transistorWidths));
    dgdwFix = dGdw - r.*(r.*dGdw);
    dwds = wMin.*(1+s./sqrt(s.^2+4))/2;
    dGds = dgdwFix.*dwds;
    
    gamma = abs((s - sPrev)'*(dGds - dGdsPrev))/(norm(dGds - dGdsPrev)^2);
    
    transistorWidthPrev = transistorWidths;
    dGdsPrev = dGds;
    sPrev = s;
    
    s = s + gamma*dGds;
    
    widBudget = @(sfix) sum(finv_s(s + sfix*ones(size(s)))) - totalWidth;
    sfix = fzero(widBudget,1);
    
    transistorWidths = finv_s(s+sfix*ones(size(s)));
    
    save(strcat('gradientState_',num2str(i)),'transistorWidths','transistorWidthPrev','s','sPrev','dGds','dGdwRaw','dGdw','gamma')

    i = i +1;
    txWidthsSave(:,i) = transistorWidths;
    fprintf('Error ||txCurrent - txPrev|| %d\n',norm(transistorWidths-transistorWidthPrev))
end

ind = nnz(saveTeola);
tCrits = saveTeola(1:ind);
txWidths = txWidthsSave(:,1:ind);

save('twoFlopSynchronizerTxWidths','tCrits','txWidths');

