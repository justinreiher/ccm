idvg = dlmread('pmos45HPtx.data');

Vd = idvg(:,1);
Vg = idvg(:,3);
Vs = idvg(:,4);
Id = idvg(:,5);


figure
hold on

for i = 1:11
    plot(Vg(1+(i-1)*101:i*101) - Vs(1+(i-1)*101:i*101),Id(1+(i-1)*101:i*101),'linewidth',1.3)
end

xlabel('V_{gs} [V]')
ylabel('I_{ds} [A]')
title('V_{gs} Characteristics Curves')
legend('V_{ds} -1.0V','V_{ds} -0.9V','V_{ds} -0.8V','V_{ds} -0.7V',...
    'V_{ds} -0.6V','V_{ds} -0.5V','V_{ds} -0.4V','V_{ds} -0.3V','V_{ds} -0.2V',...
    'V_{ds} -0.1V','V_{ds} 0.0V')

set(gca,'FontSize',18)
set(gca,'linewidth',1.2)