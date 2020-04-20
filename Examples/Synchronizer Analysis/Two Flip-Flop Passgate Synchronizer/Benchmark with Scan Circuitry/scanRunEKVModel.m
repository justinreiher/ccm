%voltage sources
gnd = vsrc('gnd',0,1);
vdd = vsrc('vdd',1,1);
clk = vtanhClock('clk',1e10,-pi/2,0.5,0.5,4,3e-10,0,1); 
clkbar = vtanhClock('clkbar',1e10,pi/2,0.5,0.5,4,3e-10,1,1);

din = vtanhDelay('din',0.5,0.5,5e10,true);

transistorWidthsScan = ones(50,1)*450e-7;

syncScan = TWO_FLOP_PG_SYNC_SCAN('syncScan',transistorWidthsScan,1);

maskScan = [syncScan.cctPath.PGFF_0.PGL_0.q,...
            syncScan.cctPath.PGFF_0.q,...
            syncScan.cctPath.PGFF_1.PGL_0.q,...
            syncScan.cctPath.q];

uVcritScan = [syncScan.cctPath.PGFF_1.PGL_0.q,...
              syncScan.cctPath.PGFF_1.PGL_0.x0];


uVcritBenchScan = zeros(1,syncScan.nodeNum);
uVcritBenchScan(uVcritScan) = [1,-1];

span = [0 1.3e-9];
tCrit = 1.2e-9;



nbScan = nestedBisection(syncScan,{syncScan.vdd,syncScan.gnd,syncScan.clk,syncScan.clkbar,syncScan.d},...
    {vdd,gnd,clk,clkbar,vdd},[5.5e-10 4e-10],{syncScan.d,din},clk,maskScan,'EKV','PTM 45nmHP');


nbAnalysisScan = nestedBisectionAnalysis(syncScan,{syncScan.vdd,syncScan.gnd,syncScan.clk,syncScan.clkbar,syncScan.d},...
    {vdd,gnd,clk,clkbar,vdd},{syncScan.d,din},clk,syncScan.cctPath.PGFF_1.PGL_0.q,...
    'EKV','PTM 45nmHP');


benchScanRuns = nbScan.simulate('2flopSyncEKVScan',span,tCrit);

[tbenchScan,gdotgBenchScan] = nbAnalysisScan.bisectionAnalysis('benchScanAnalysisEKV',benchScanRuns,uVcritBenchScan,tCrit);



figure
hold on
yyaxis left
plot(tbenchScan,gdotgBenchScan,'linewidth',1.2)

ylabel('Gain [1/s]')
xlabel('Time [s]')
axis([6.5e-10 9e-10 -0.5e11 2e11])
yyaxis right
plot(tbenchScan,clk.V(tbenchScan),'-g','linewidth',1.2)
ylabel('Voltage [V]')
set(gca,'FontSize',14)
legend('Scan','Clk')
title('Instantaneous Gain $\lambda(t)$ EKV model','interpreter','latex','FontSize',20)
set(gca,'fontsize',18)