% Fitting script, 
% Yan Peng, UBC, 2016/10/19
% Justin Reiher, 2017/05/30 modifications

clear all

ids_sweepVds_N = './45nmHP/nmos45HP_VdsSweep.txt';
ids_sweepVgs_N = './45nmHP/nmos45HP_VgsSweep.txt';
ids_sweepVgd_N = './45nmHP/nmos45HP_VgdSweep.txt';
ids_sweepVsd_N = './45nmHP/nmos45HP_VsdSweep.txt';

ids_sweepVds_P = './45nmHP/pmos45HP_VdsSweep.txt';
ids_sweepVgs_P = './45nmHP/pmos45HP_VgsSweep.txt';
ids_sweepVgd_P=  './45nmHP/pmos45HP_VgdSweep.txt';
ids_sweepVsd_P = './45nmHP/pmos45HP_VsdSweep.txt';


% read in the data for Ids vs Vds and Ids vs Vgs
%Constants for Ids vs Vds curves
Vx_start =0;
Vx_step_size = 0.01;
Vx_stop = 1.0;

Vy_start = 0;
Vy_step_size = 0.1;
Vy_stop = 1.0;

Vmax = 1.0;

sweeps = [Vmax,Vx_start,Vx_step_size,Vx_stop,Vy_start,Vy_step_size,Vy_stop];
	
	
%% mos transistor parameters
mosParams.W = 450E-7;  	%Transistor Width in [cm]
mosParams.Lg = 45e-7;	 	%Transistor Gate Lenght in [cm]
mosParams.alpha = 13; 		%Scaling Term
mosParams.beta = 0.15;         %Drain Induced Barrier Lowering term 
mosParams.Vth = 0.5;		%The threshold Voltage
mosParams.Io = 1; 
mosParams.gamma = 0.2;    %Paramaters to model body effect
mosParams.phi = 0.5;
mosParams.type = 1; %1 denotes an nmos device, -1 is a pmos device
%prompt = 'Maybe you want to look at the plots first? (true or false)\n';
%x = input(prompt);
%if isempty(x)
%    to_plot = false;
%else
%    to_plot = true;
%end

[idvd_n, idvg_n, idvs_n, idvd_p, idvg_p, idvs_p] = read_data(ids_sweepVds_N, ids_sweepVsd_N, ids_sweepVg_N, ids_sweepVs_N,...
                                                            ids_sweepVd_P,ids_sweepVg_P,ids_sweepVs_P, sweeps);


mosParams.type = 1;
mosParamsN = fit_ekv(mosParams,idvd_n,idvg_n,idvs_n,sweeps,'n');

mosParams.type = -1;
mosParamsP = fit_ekv(mosParams,idvd_p,idvg_p,idvs_p,sweeps,'p');

plot_ekv(mosParamsN,sweeps,idvd_n,idvg_n,idvs_n,false,'n')
plot_ekv(mosParamsP,sweeps,idvd_p,idvg_p,idvs_p,false,'p')
