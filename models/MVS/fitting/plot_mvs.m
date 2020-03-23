function []=plot_mvs(input_parms, coeff_op_final, Vy, IdVd, IdVg, take_log)
%% plotting output data
%%% Read output curve data for plotting
%Vy is an array of scaling data in the form of:
%
% Vy = [Vds_start, Vds_step_size, Vds_stop, Vgs_step_start, Vgs_steps, Vgs_step_stop,...
% 	Vds_step_start, Vds_step_stop, Vgs_start, Vgs_step_size, Vgs_stop];
clear Id_data
clear bias_data
clear Vg* Id_out

%%% ####### 
Vystart = Vy(1);
Vystep = Vy(2);
Vyend = Vy(3);
Vypre = Vystart:Vystep:Vyend;
vv=length(Vypre);
Vy_data=IdVd(:,1);
Vg_data=IdVd(:,2);
Id_data = IdVd(:,end);
bias_data(:,1)=Vy_data;
bias_data(:,2)=Vg_data;
Vb=0;
bias_data(:,3)=Vb;
Vx=0;
bias_data(:,4)=Vx;

[Id_optim1] = mvs_si_1_1_0(input_parms, coeff_op_final, bias_data);
ll = length(Id_optim1)/vv;
figure
hold on
legendInfo = cell(1,ll*2);
for len_vd = 1:ll
 color = rand(1,3);
 plot(Vypre,(Id_data((len_vd-1)*vv+1:len_vd*vv)),'o','color',color);
 plot(Vypre,10.^(Id_optim1((len_vd-1)*vv+1:len_vd*vv)),'color',color);
 hold on
legendInfo{2*len_vd-1} = ['Vgs = ' num2str(Vystep*len_vd)];
 legendInfo{2*len_vd} = ['Vgs = ' num2str(Vystep*len_vd)];
end
%columnlegend(6,legendInfo,'NorthWest');
xlabel('Vds');
ylabel('Ids');
title('Sweeping Vds for each Vgs');


%% plotting transfer
%% Read transfer curve data for plotting 
clear bias_data
clear Id_data

Vgpre=IdVg(:,1);
IdVg_loVd=IdVg(:,2);
IdVg_hiVd=IdVg(:,3);
Vy_data=[Vystart*ones(1,length(Vgpre))';Vyend*ones(1,length(Vgpre))'];
Vg_data=[Vgpre;Vgpre];
bias_data(:,1)=Vy_data;
bias_data(:,2)=Vg_data;
Vb=0;
bias_data(:,3)=Vb;
Vx=0;
bias_data(:,4)=Vx;
Id_optim2=mvs_si_1_1_0(input_parms, coeff_op_final, bias_data);
figure
if take_log
    semilogy(Vgpre,IdVg_loVd,'bo',Vgpre,IdVg_hiVd,'ro');
else
    plot(Vgpre,IdVg_loVd,'bo',Vgpre,IdVg_hiVd,'ro');
end
hold on
if take_log
    plot(Vgpre,10.^(Id_optim2(1:length(Vgpre))),'b',Vgpre,10.^(Id_optim2(length(Vgpre)+1:end)),'r');
	title('Log plot of sweeping Vgs for two Vds');
else
    plot(Vgpre,10.^(Id_optim2(1:length(Vgpre))),'b',Vgpre,10.^(Id_optim2(length(Vgpre)+1:end)),'r');
	title('Plot of sweeping Vgs for two Vds');
end
legend(['Vds = ' num2str(Vystart) 'V'], ['Vds = ' num2str(Vyend) 'V']);
xlabel('Vgs');
ylabel('Ids');


end
