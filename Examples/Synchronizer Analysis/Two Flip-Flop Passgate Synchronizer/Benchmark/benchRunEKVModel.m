%voltage sources
gnd = vsrc('gnd',0,1);
vdd = vsrc('vdd',1,1);
clk = vtanhClock('clk',1e10,-pi/2,0.5,0.5,4,3e-10,0,1); 
clkbar = vtanhClock('clkbar',1e10,pi/2,0.5,0.5,4,3e-10,1,1);

din = vtanhDelay('din',0.5,0.5,5e10,true);

transistorNumbering = 1:1:42;
nDevices = odd(transistorNumbering);
pDevices = even(transistorNumbering);

transistorWidths = ones(42,1)*450e-7;

synch = TWO_FLOP_PG_SYNC('sync',transistorWidths,1);

mask = [synch.cctPath.PGFF_0.PGL_0.q,...
        synch.cctPath.PGFF_0.q,...
        synch.cctPath.PGFF_1.PGL_0.q,...
        synch.cctPath.q];

uVcrit = [synch.cctPath.PGFF_1.PGL_0.q,...
          synch.cctPath.PGFF_1.PGL_0.x0];

uVcritBench = zeros(1,synch.nodeNum);
uVcritBench(uVcrit) = [1,-1];

span = [0 1.3e-9];
tCrit = 1.2e-9;

nbBench = nestedBisection(synch,{synch.vdd,synch.gnd,synch.clk,synch.clkbar,synch.d},...
    {vdd,gnd,clk,clkbar,vdd},[5.5e-10 4e-10],{synch.d,din},clk,mask,'EKV','PTM 45nmHP');

nbAnalysis = nestedBisectionAnalysis(synch,{synch.vdd,synch.gnd,synch.clk,synch.clkbar,synch.d},...
    {vdd,gnd,clk,clkbar,vdd},{synch.d,din},clk,synch.cctPath.PGFF_1.PGL_0.q,'EKV','PTM 45nmHP');


benchRuns = nbBench.simulate('2flopSyncEKV',span,tCrit);

[tbench,gdotgBench] = nbAnalysis.bisectionAnalysis('benchAnalysisEKV',benchRuns,uVcritBench,tCrit);


figure
hold on
yyaxis left
plot(tbench,gdotgBench,'-m','linewidth',1.2)

ylabel('\lambda(t) [1/s]')
xlabel('Time [s]')
axis([dataTransition tbench(end) -inf inf])
yyaxis right
plot(tbench,clk.V(tbench),'-g','linewidth',1.2)
ylabel('Voltage [V]')
set(gca,'FontSize',14)
legend('\lambda(t) - Benchmark','Clock')
title('Instantaneous Gain $\lambda(t)$ EKV model','interpreter','latex','FontSize',20)
set(gca,'fontsize',18)

load('benchAnalysisEKV')

gBench = zeros(1,length(tbench));
for i = 1:length(tbench)
    gBench(i) = wNorm_t(tbench(i))'*beta_t(tbench(i));
end