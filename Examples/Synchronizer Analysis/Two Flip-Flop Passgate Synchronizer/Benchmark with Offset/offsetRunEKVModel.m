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


transistorWidthsOffset = transistorWidths;
transistorWidthsOffset(8) = 600e-7;
transistorWidthsOffset(7) = 300e-7;


syncOff = TWO_FLOP_PG_SYNC('sync',transistorWidthsOffset,1);

mask = [syncOff.cctPath.PGFF_0.PGL_0.q,...
        syncOff.cctPath.PGFF_0.q,...
        syncOff.cctPath.PGFF_1.PGL_0.q,...
        syncOff.cctPath.q];

uVcrit = [syncOff.cctPath.PGFF_1.PGL_0.q,...
          syncOff.cctPath.PGFF_1.PGL_0.x0];

uVcritOffset = zeros(1,syncOff.nodeNum);
uVcritOffset(uVcrit) = [1,-1];

span = [0 1.3e-9];
tCrit = 1.2e-9;


nbOffset = nestedBisection(syncOff,{syncOff.vdd,syncOff.gnd,syncOff.clk,syncOff.clkbar,syncOff.d},...
    {vdd,gnd,clk,clkbar,vdd},[5.5e-10 4e-10],{syncOff.d,din},clk,mask,'EKV','PTM 45nmHP');


nbAnalysisOff = nestedBisectionAnalysis(syncOff,{syncOff.vdd,syncOff.gnd,syncOff.clk,syncOff.clkbar,syncOff.d},...
   {vdd,gnd,clk,clkbar,vdd},{syncOff.d,din},clk,syncOff.cctPath.PGFF_1.PGL_0.q,'EKV','PTM 45nmHP');

offsetRuns = nbOffset.simulate('2flopSyncOffsetEKV',span,tCrit);

[toffset,gdotgOffset] = nbAnalysisOff.bisectionAnalysis('offsetAnalysisEKV',offsetRuns,uVcritOffset,tCrit);


figure
hold on
yyaxis left
plot(toffset,gdotgOffset,'-m','linewidth',1.2)

ylabel('\lambda(t) [1/s]')
xlabel('Time [s]')
axis([dataTransition toffset(end) -inf inf])
yyaxis right
plot(toffset,clk.V(toffset),'-g','linewidth',1.2)
ylabel('Voltage [V]')
set(gca,'FontSize',14)
legend('\lambda(t) - Offset','Clock')
title('Instantaneous Gain $\lambda(t)$ EKV model','interpreter','latex','FontSize',20)
set(gca,'fontsize',18)

load('offsetAnalysisEKV')

gOffset = zeros(1,length(toffset));
for i = 1:length(toffset)
    gOffset(i) = wNorm_t(toffset(i))'*beta_t(toffset(i));
end