function [coeff_op_final,idvd,idvg,input_params] = fit_mvs(mos_params,idvd,idvg,outfile_name)
% Rewritten by Yan Peng, UBC, 2016/10/19

% choose_mos = 1;             % choose_mos = 1 for nmos and 0 for pmos
Tjun=300;                   % Temperature in degrees Kelvin

 type = mos_params(8);    % *** nFET: type=1  pFET: type=-1

% 16nm mos
W= mos_params(1)          % Note: 1e-7cm = 1nm
Lgdr = mos_params(2);
dLg = mos_params(3); 
Cg = mos_params(4);

%% UNCOMMENT following 4 lines for 32-nm data
% W=1e-4;                     % *** Width [cm]
% Lgdr = 32e-7;                 % *** Gate length [cm]
% dLg= 9e-7;                  % *** dLg=L_g-L_c *1e-7 (default 0.3xLg_nom)
% Cg=2.57e-6;                 %*** Gate capacitance in F/cm2
 
%% UNCOMMENT following 4 lines for 45-nm data
% W=1e-4;                     % *** Width [cm] 
% Lgdr = 45e-7;               % *** Gate length [cm] 
% dLg= 7.56e-7;               % *** dLg=L_g-L_c*1e-7 (default 0.3xLg_nom) 
% Cg=2.55e-6;                 % *** Gate capacitance in F/cm^2
Vy_max = mos_params(5);
beta = mos_params(7);
alpha=3.5;                  % *** Charge Vt trasistion factor (don't change this mostly between 3.0 and 4.0)

% following parameters are not used for optimization routine at all. But
% they must be provided since they are used as inputs.
IdA=1.5E-9;               % *** Desired Id at Vg=VgA and Vd=VdA: Must be well in sub-threshold
VdA=0.05;                   % *** Vd [V] corresponding to IdA
VgA=0.1;                    % *** Vg [V] corresponding to IdA 
                            % *** Above values overriden if Vt0 is
                            % specified directly

etov = 1.75e-007;          % *** Equivalent thickness of dielectric at S/D-G overlap [cm]
Cif = mos_params(6);              % *** Inner fringing S or D capacitance [F/cm] 
Cof = mos_params(6);              % *** Outer fringing S or D capacitance [F/cm]

rv=1.0;                     % *** Ratio vxo(strong inversion)/vxo(weak inversion)
                            % *** Set rv=1 for constant vxo (zeta
                            % irrelevant but do not set zeta=0)
zeta=1;                     % *** Parameter determines transtion Vg for vxo 


phib=1.2;                   % *** ~abs(2*phif) [V]
gamma=0.2;                  % *** Body factor  [sqrt(V)]
mc=0.2;                     % *** Carrier effective mass, raltive to m_0. 
                            % *** Choose an appropriate value between 0.01
                            % *** to 10. For, values outside of this range,
                            % *** convergence or accuracy of results is not
                            % *** guaranteed.
                            
CTM_select = 1;             % *** Parameter to select charge-transport model 
                            % *** if CTM_select = 1, then classic DD-NVSAT
                            % *** model is used; for CTM_select other than
                            % *** 1,blended DD-NVSAT and ballistic charge
                            % *** transport model is used.

%%% define global inputs
input_params =[type;W;Lgdr;dLg;gamma;phib;Cg;Cif;Cof;etov;mc;Tjun;beta;alpha;CTM_select;zeta]; % set of input parameters 

%%%% *********

% intial guess for parameters (a total of 7 parameters that need to be
% fitted since we assume Rs0 = Rd0.) the condition Rs0 = Rd0 has been
% hard-coded in the model file. If Rs0 and Rd0 are different, then the
% model file must be appropritaely tweaked.

Rs0=100;                     % *** Access resistance for terminal "x" [ohm-micron] (Typically Rs)  
Rd0=100;                     % *** Access resistance for terminal "y" (Typically assume Rs=Rd)


delta=0.15;                 % *** DIBL [V/V] 
S0=0.1;                     % *** Subthreshold swing at T=27 C and Vd=VdA [V/decade]
phit = 8.617e-5*(273+Tjun);    % kT/q
S=S0*(Tjun+273)/300;
n = S/(log(10)*phit);
nd=0;                       % *** Factor allowing for modest punchthrough.  
                            % *** Normally, nd=0.  If some punchtrhough
                            % 0<nd<0.4
n0=n-nd*VdA;                % Intrinsic swing n-factor at Tjun

vxo=1.2;                    % *** Virtual source velocity [cm/s]    
mu=200;                     % *** Mobility [cm^2/V.s]

%Vt0 =
%VT_new_SR(W,Lg,dLg,IdA,VgA,VdA,Cg,delta,n0,nd,vxo,rv,zeta,mu,phit,alpha,beta);
%%% initial computation of Vt0 (overwritten)
Vt0 = 0.4; % hardcoded Vt0. will be optimized in the process


coeff_init=[Rs0; Rd0; delta; n0; nd; vxo; mu; Vt0]; % matrix of coefficients to be optimized. 7 coefficients. Rd0 is dummy.

[coeff_op_tran] = optimize_transfer(input_params,Vy_max,coeff_init, idvg); % first iteration
[coeff_op_out] = optimize_output(input_params,coeff_op_tran, idvd); % first iteration


iter=1; % this can be changed to a suitable iteration. optimized params from transfer are fed as initial guesses into output chars.
for len_iter=2:iter
[coeff_op_out(:,len_iter)] = optimize_output(input_params,coeff_op_tran(:,len_iter-1), idvg);
[coeff_op_tran(:,len_iter)] = optimize_transfer(input_params,Vy_max,coeff_op_out(:,len_iter-1), idvd);
end
coeff_op_avg = (coeff_op_out+coeff_op_tran)/2; % there may be a bit of discrepancy between the transfer and output and therfore, final optimized params are taken as average.
coeff_op_final = mean(coeff_op_avg(:,max(1,end-8):end),2);
% coeff_op_final = [1.1; 1.1; ones(length(coeff_op_final)-2,1)] .* coeff_op_final
%%% ########################################################

% writing optimized coefficients in the output_text.txt file.
format long e
row_dim = length(coeff_init);
col_dim = iter;
format_string = cell(col_dim,1);
for nC = 1:col_dim
    format_string{nC} = ['%12.4f'];
end

outfile_name = [outfile_name, '.coeff'];
fileID = fopen(outfile_name, 'w');
fprintf(fileID, 'coeff_op_tran');
fprintf(fileID, '\n');

%Coeff_op_tran
row_labels_tran = cell(row_dim,1);
row_labels_tran{1} = ['Rs0'];
row_labels_tran{2} = ['Rd0'];        % this parameter does not matter.
row_labels_tran{3} = ['DIBL'];
row_labels_tran{4} = ['n0'];
row_labels_tran{5} = ['nd'];
row_labels_tran{6} = ['vxo'];
row_labels_tran{7} = ['mu'];
row_labels_tran{8} = ['Vt0'];

for nA = 1:row_dim
    fprintf(fileID, '%8s\t', row_labels_tran{nA});
    fprintf(fileID, [format_string{:}, '\n'], coeff_op_tran(nA,:));
end
fprintf(fileID, '\n');

%Coeff_op_out
row_labels_out = cell(row_dim,1);
row_labels_out{1} = ['Rs0'];
row_labels_out{2} = ['Rd0'];  % this parameter does not matter.
row_labels_out{3} = ['DIBL'];
row_labels_out{4} = ['n0'];
row_labels_out{5} = ['nd'];
row_labels_out{6} = ['vxo'];
row_labels_out{7} = ['mu'];
row_labels_out{8} = ['Vt0'];

fprintf(fileID, 'coeff_op_out');
fprintf(fileID, '\n');
for nB = 1:row_dim
    fprintf(fileID, '%8s\t', row_labels_out{nB});
    fprintf(fileID, [format_string{:}, '\n'], coeff_op_out(nB,:));

end

%% average optimized coefficients
fprintf(fileID, '\n');
row_labels_avg = cell(row_dim,1);
row_labels_avg{1} = ['Rs0'];
row_labels_avg{2} = ['Rd0'];     % this parameter does not matter
row_labels_avg{3} = ['DIBL'];
row_labels_avg{4} = ['n0'];
row_labels_avg{5} = ['nd'];
row_labels_avg{6} = ['vxo'];
row_labels_avg{7} = ['mu'];
row_labels_avg{8} = ['Vt0'];

fprintf(fileID, 'coeff_op_avg');
fprintf(fileID, '\n');
for navg=1:row_dim
    fprintf(fileID, '%8s\t', row_labels_avg{navg});
    fprintf(fileID, [format_string{:}, '\n'], coeff_op_avg(navg,:));
end


fprintf(fileID, '\n');
row_labels_avg = cell(row_dim,1);
row_labels_avg{1} = ['Rs0'];
row_labels_avg{2} = ['Rd0'];     % this parameter does not matter
row_labels_avg{3} = ['DIBL'];
row_labels_avg{4} = ['n0'];
row_labels_avg{5} = ['nd'];
row_labels_avg{6} = ['vxo'];
row_labels_avg{7} = ['mu'];
row_labels_avg{8} = ['Vt0'];

fprintf(fileID, 'coeff_op_final');
fprintf(fileID, '\n');
for nfinal=1:row_dim
    fprintf(fileID, '%8s\t', row_labels_avg{nfinal});
    fprintf(fileID, [format_string{:}, '\n'], coeff_op_final(nfinal,:));
end

fprintf(fileID, '\n');
fclose(fileID);
%%% ########################################################
%%% end of file that contains data.

end
