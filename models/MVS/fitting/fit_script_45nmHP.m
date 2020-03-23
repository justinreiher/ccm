% Fitting script, 
% Yan Peng, UBC, 2016/10/19
% Justin Reiher, 2017/05/30 modifications

clear all

idvd_n_file = './45nmHP/nmos45nmHP_IV.data';
idvg_n_file = './45nmHP/nmos45nmHP_transfer.data';

idvd_p_file = './45nmHP/pmos45nmHP_IV_2.data';
idvg_p_file = './45nmHP/pmos45nmHP_transfer_2.data';

outfilename_n = './45nmHP_n';
outfilename_p = './45nmHP_p_2';

% read in the data for Ids vs Vds and Ids vs Vgs
%Constants for Ids vs Vds curves
Vds_start =0.05;
Vds_step_size = 0.01;
Vds_stop = 1.0;

Vgs_step_start = 0.3;
Vgs_steps = 0.1;
Vgs_step_stop = 1.0;

%Constants for Ids vs Vgs transfer curves

Vds_step_start = 0.05;
Vds_step_stop = 1.0;

Vgs_start = 0;
Vgs_step_size =0.01;
Vgs_stop =1.0;

% [ Vds_start, Vds_step_size, Vds_stop, Vgs_step_start, Vgs_steps, step_Vgs_stop,Vds_step_start, Vds_step_end Vgs_start, Vgs_step, Vgs_stop]

scale = [Vds_start, Vds_step_size, Vds_stop, Vgs_step_start, Vgs_steps, Vgs_step_stop,...
 	Vds_step_start, Vds_step_stop, Vgs_start, Vgs_step_size, Vgs_stop];
	
	
%% mos transistor parameters
mos_params = zeros(1,5);
mos_params(1) = 1000E-7;  		%Transistor Width in [cm]
mos_params(2) = 45e-7;	 		%Transistor Gate Lenght in [cm]
mos_params(3) = 3.75e-7; 		%dLg = Lg - Lc * 1e-7 (default is 0.3*LgNom) in [cm]
mos_params(4) = 2.837e-6; 	%Gate capacitance in [F/cm^2]
mos_params(5) = Vgs_stop;		%The upper bound on the voltage sweep
mos_params(6) = 1.9856e-12; 
%prompt = 'Maybe you want to look at the plots first? (true or false)\n';
%x = input(prompt);
%if isempty(x)
%    to_plot = false;
%else
%    to_plot = true;
%end

[idvd_n, idvd_p, idvg_n, idvg_p] = read_data(idvd_n_file, idvd_p_file, idvg_n_file, idvg_p_file, scale);


% Fit the data
coeff = zeros(8,2);
mos_params(7) = 1.8; % *** Saturation factor. Typ. nFET=1.8, pFET=1.6
mos_params(8)=1;
[coeff(:,1), IdVd_n, IdVg_n, input_parms_n] = fit_mvs(mos_params, idvd_n, idvg_n,outfilename_n);
plot_mvs(input_parms_n, coeff(:,1), scale, IdVd_n, IdVg_n, false);
drawnow
mos_params(7) = 1.6; % *** Saturation factor. Typ. nFET=1.8, pFET=1.6
%mos_params(7) = -1; %fitting with this flag set does not work
[coeff(:,2), IdVd_p, IdVg_p, input_parms_p] = fit_mvs(mos_params, abs(idvd_p), abs(idvg_p),outfilename_p);
plot_mvs(input_parms_p, coeff(:,2), scale, IdVd_p, IdVg_p, false);
drawnow

