filename = '3stageRingOsc45nm.data';

Vdd = vsrc('vdd',1.0,1);
Gnd = vsrc('gnd',0.0,1);

oscData = dlmread(filename);
timeSim = oscData(:,1)-200e-12;
Vo3 = oscData(:,2);
Vo2 = oscData(:,4);
Vo1 = oscData(:,6);

ringOsc = THREE_RING_OSC('osc',450e-7,1);
myDict = modelDictionary;
tbOptionsFull = struct('capModel','full','capScale',1,'vdd',1,'temp',298,'numParallelCCTs',1,'debug',false);
tbOptionsGnd  = struct('capModel','gnd','capScale',1,'vdd',1,'temp',298,'numParallelCCTs',1,'debug',false);


tb_ekvFull = testbench(ringOsc,{ringOsc.vdd,ringOsc.gnd},{Vdd,Gnd},myDict,'EKV','PTM 45nmHP',tbOptionsFull);
tb_ekvGnd  = testbench(ringOsc,{ringOsc.vdd,ringOsc.gnd},{Vdd,Gnd},myDict,'EKV','PTM 45nmHP',tbOptionsGnd);

tb_mvsFull  = testbench(ringOsc,{ringOsc.vdd,ringOsc.gnd},{Vdd,Gnd},myDict,'MVS','PTM 45nmHP',tbOptionsFull);
tb_mvsGnd  = testbench(ringOsc,{ringOsc.vdd,ringOsc.gnd},{Vdd,Gnd},myDict,'MVS','PTM 45nmHP',tbOptionsGnd);

[tEkvFull,VekvFull] = tb_ekvFull.simulate([0 2e-10],[0,0,0,1,0]);
[tEkvGnd,VekvGnd]   = tb_ekvGnd.simulate([0 2e-10],[0,0,0,1,0]);
[tMvsFull,VmvsFull]   = tb_mvsFull.simulate([0 2e-10],[0,0,0,1,0]);
[tMvsGnd,VmvsGnd]   = tb_mvsGnd.simulate([0 2e-10],[0,0,0,1,0]);



figure
subplot(4,1,1)
hold on
plot(timeSim,Vo3,'o')
plot(tEkvFull-1.8e-11-1.1e-12,VekvFull(:,5))
plot(tEkvFull,0.5*ones(1,length(tEkvFull)),'r--')
title('Three Ring Oscillator Simulation')
ylabel('Output Voltage [V]')
legend('ngSpice','EKV Model - Full')
set(gca,'fontsize',18)
axis([0,timeSim(end),-inf,inf])
set(findall(gca,'type','line'),'linewidth',1.2)

subplot(4,1,2)
hold on
plot(timeSim,Vo3,'o')
plot(tEkvGnd-1.8e-11-1.1e-12,VekvGnd(:,5))
plot(tEkvGnd,0.5*ones(1,length(tEkvGnd)),'r--')
ylabel('Output Voltage [V]')
legend('ngSpice','EKV Model - Gnd')
set(gca,'fontsize',18)
axis([0,timeSim(end),-inf,inf])
set(findall(gca,'type','line'),'linewidth',1.2)

subplot(4,1,3)
hold on
plot(timeSim,Vo3,'o')
plot(tMvsFull-1.8e-11-1.1e-12,VmvsFull(:,5))
plot(tMvsFull,0.5*ones(1,length(tMvsFull)),'r--')
ylabel('Output Voltage [V]')
legend('ngSpice','MVS Model - Full')
set(gca,'fontsize',18)
axis([0,timeSim(end),-inf,inf])
set(findall(gca,'type','line'),'linewidth',1.2)

subplot(4,1,4)
hold on
plot(timeSim,Vo3,'o')
plot(tMvsGnd-1.8e-11-1.1e-12,VmvsGnd(:,5))
plot(tMvsGnd,0.5*ones(1,length(tMvsGnd)),'r--')
xlabel('Time [s]')
ylabel('Output Voltage [V]')
legend('ngSpice','MVS Model - Gnd')
set(gca,'fontsize',18)
axis([0,timeSim(end),-inf,inf])
set(findall(gca,'type','line'),'linewidth',1.2)