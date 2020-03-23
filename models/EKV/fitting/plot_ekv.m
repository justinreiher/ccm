function plot_ekv(input_parms, sweeps, IdVd, IdVg, IdVs, take_log,type)
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
Vystart = sweeps(2);
Vystep = sweeps(3);
Vyend = sweeps(4);
Vypre = Vystart:Vystep:Vyend;
vv = length(Vypre);
Vy1 = IdVd(:,1)';
Vy2 = IdVs(:,1)';

Vg1 = IdVd(:,2)';
Vg2 = IdVs(:,2)';
 
Vx1 = IdVd(:,3)';
Vx2 = IdVs(:,3)';

Vb1 = IdVd(:,4)';
Vb2 = IdVs(:,4)';

Id_data = IdVd(:,5)';
Id_data2 = IdVs(:,5)';

IdsFit1 = ekv_Ids(input_parms, [Vy1;Vg1;Vx1;Vb1]);
IdsFit2 = ekv_Ids(input_parms, [Vy2;Vg2;Vx2;Vb2]);
runs = length(IdsFit1)/vv;
figure
hold on
legendInfo = cell(1,11*2);
for len_vd = 1:runs
 color = rand(1,3);
 start = (len_vd-1)*vv+1;
 stop  = len_vd*vv;
 Ids1 = Id_data(start:stop);
 Ids2 = Id_data2(start:stop);
 
 Vds1 = IdVd(start:stop,1)-IdVd(start:stop,3);
 Vds2 = IdVs(start:stop,1)-IdVs(start:stop,3);
 
 plot(Vds1,Ids1,'o','color',color);
 plot(Vds2,Ids2,'o','color',color);
 plot(Vds1,IdsFit1(start:stop),'color',color);
 plot(Vds2,IdsFit2(start:stop),'color',color);
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

mid = round(length(IdVg)/(2*vv));
endPoint = length(IdVg)/vv;

lowStart = 1+vv;
lowEnd   = 2*vv;

midStart = 1+(mid-1)*vv;
midEnd   = mid*vv;

highStart     = 1+(endPoint -2)*vv;
highEnd  = (endPoint-1)*vv;


IdVg_loVd=IdVg(lowStart:lowEnd,5);
IdVg_midVd = IdVg(midStart:midEnd,5);
IdVg_hiVd=IdVg(highStart:highEnd,5);

Vy = IdVg(:,1)';
Vg = IdVg(:,2)';
Vx = IdVg(:,3)';
Vb = IdVg(:,4)';


Id_optim2=ekv_Ids(input_parms, [Vy;Vg;Vx;Vb]);
figure
if take_log
    semilogy(Vg(lowStart:lowEnd)-Vx(lowStart:lowEnd),input_parms.type*IdVg_loVd,'bo',...
        Vg(midStart:midEnd)-Vx(midStart:midEnd),input_parms.type*IdVg_midVd,'go',...
        Vg(highStart:highEnd)-Vx(highStart:highEnd),input_parms.type*IdVg_hiVd,'ro');
else
    plot(Vg(lowStart:lowEnd)-Vx(lowStart:lowEnd),input_parms.type*IdVg_loVd,'bo',...
        Vg(midStart:midEnd)-Vx(midStart:midEnd),input_parms.type*IdVg_midVd,'go',...
        Vg(highStart:highEnd)-Vx(highStart:highEnd),input_parms.type*IdVg_hiVd,'ro');
end
hold on
if take_log
    semilogy(Vg(lowStart:lowEnd)-Vx(lowStart:lowEnd),input_parms.type*Id_optim2(lowStart:lowEnd),'b',...
        Vg(midStart:midEnd)-Vx(midStart:midEnd),input_parms.type*Id_optim2(midStart:midEnd),'g',...
        Vg(highStart:highEnd)-Vx(highStart:highEnd),input_parms.type*Id_optim2(highStart:highEnd),'r');
	title('Log plot of sweeping Vgs for two Vds');
else
    plot(Vg(lowStart:lowEnd)-Vx(lowStart:lowEnd),input_parms.type*Id_optim2(lowStart:lowEnd),'b',...
        Vg(midStart:midEnd)-Vx(midStart:midEnd),input_parms.type*Id_optim2(midStart:midEnd),'g',...
        Vg(highStart:highEnd)-Vx(highStart:highEnd),input_parms.type*Id_optim2(highStart:highEnd),'r');
	title('Plot of sweeping Vgs for two Vds');
end
lo = Vy(lowStart:lowEnd) - Vx(lowStart:lowEnd);
mid = Vy(midStart:midEnd) - Vx(midStart:midEnd);
high = Vy(highStart:highEnd) - Vx(highStart:highEnd);
legend(['Vds = ' num2str(lo(1)) 'V'], ['Vds = ' num2str(mid(1)) 'V'], ['Vds = ' num2str(high(1)) 'V']);
xlabel('Vgs');
ylabel('Ids');


end
