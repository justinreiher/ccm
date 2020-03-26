%voltage sources
gnd = vsrc('gnd',0,1);
vdd = vsrc('vdd',1,1);
clk = vtanhClock('clk',1e10,-pi/2,0.5,0.5,4,3e-10,0,1); 
clkbar = vtanhClock('clkbar',1e10,pi/2,0.5,0.5,4,3e-10,1,1);

din = vtanhDelay('din',0.5,0.5,5e10,true);

myDict = modelDictionary;
transistorNumbering = 1:1:12;
nDevices = odd(transistorNumbering);
pDevices = even(transistorNumbering);

transistorWidths = ones(12,1)*450e-7;

pgREF = PASSGATE_LATCH_TEST('pgL',450e-7,1);
pg = PASSGATE_LATCH_TEST('pgL',transistorWidths,1);
mask = [6,6];

uVcrit = [6,9];

uVcritBench = zeros(1,10);
uVcritBench(uVcrit) = [1,-1];

span = [0 7e-10];
tCrit = 6.5e-10;

nbPGLatch = nestedBisection(pg,{pg.vdd,pg.gnd,pg.clk,pg.clkbar,pg.d},...
    {vdd,gnd,clk,clkbar,vdd},[5.5e-10 4e-10],{pg.d,din},clk,mask,uVcrit,myDict,'EKV','PTM 45nmHP');

nbPGLAnalysis = nestedBisectionAnalysis(pg,{pg.vdd,pg.gnd,pg.clk,pg.clkbar,pg.d},...
    {vdd,gnd,clk,clkbar,vdd},{pg.d,din},clk,mask(end),myDict,'EKV','PTM 45nmHP');


benchRuns = nbPGLatch.simulate('pgLatchEKV',span,tCrit);

[tpgLatch,lambdaPGL] = nbPGLAnalysis.bisectionAnalysis('pgLatchAnalysisEKV',benchRuns,uVcritBench,tCrit);

midClk = @(t) clk.V(t) - 0.5;

tOff = fzero(midClk,5e-10);

load('pgLatchAnalysisEKV')

figure
hold on
yyaxis left
plot(tpgLatch-tOff,lambdaPGL,'-m','linewidth',1.2)

ylabel('\lambda(t) [1/s]')
xlabel('Time [s]')
axis([dataTransition-tOff tpgLatch(end)-tOff -inf inf])
yyaxis right
plot(tpgLatch-tOff,clk.V(tpgLatch),'-g','linewidth',1.2)
ylabel('Voltage [V]')
set(gca,'FontSize',14)
legend('\lambda(t) - PG Latch','Clock')
title('Instantaneous Gain $\lambda(t)$ MVS model','interpreter','latex','FontSize',20)
set(gca,'fontsize',18)


gPgLatch = zeros(1,length(tpgLatch));
for i = 1:length(tpgLatch)
    gPgLatch(i) = wNorm_t(tpgLatch(i))'*beta_t(tpgLatch(i));
end