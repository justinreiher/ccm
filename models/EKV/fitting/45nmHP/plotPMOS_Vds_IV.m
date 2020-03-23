idvd = dlmread('pmos45HP.data');
idvg = dlmread('pmos45HPtx.data');

Vd = idvd(:,1);
Vg = idvd(:,3);
Vs = idvd(:,4);
Id = idvd(:,5);


figure
hold on

for i = 1:11
    plot(Vd(1+(i-1)*101:i*101) - Vs(1+(i-1)*101:i*101),Id(1+(i-1)*101:i*101),'linewidth',1.3)
end

xlabel('V_{ds} [V]')
ylabel('I_{ds} [A]')
title('V_{ds} Characteristics Curves')
legend('V_{gs} -1.0V','V_{gs} -0.9V','V_{gs} -0.8V','V_{gs} -0.7V',...
    'V_{gs} -0.6V','V_{gs} -0.5V','V_{gs} -0.4V','V_{gs} -0.3V','V_{gs} -0.2V',...
    'V_{gs} -0.1V','V_{gs} 0.0V')

set(gca,'FontSize',18)
set(gca,'linewidth',1.2)