%% optimization script for transfer data set.
% written on July 08, 2013
% Author: Shaloo Rakheja, MIT
% Rewritten by Yan Peng, UBC, 2016/10/19


function [coeff_op_tran] = optimize_transfer(input_params,Vmax, coeff_init, IdVg)
clear Id_data
clear bias_data
clear Vg* Vy* 

IdVg = abs(IdVg);

Vgpre=IdVg(1:end,1);
IdVg_loVd=IdVg(1:end,2);
IdVg_hiVd=IdVg(1:end,3);
Vy_data=[0.05*ones(1,length(Vgpre))';Vmax*ones(1,length(Vgpre))'];
Vg_data=[Vgpre;Vgpre];
Id_data(:,1)=[IdVg_loVd;IdVg_hiVd];
bias_data(:,1)=Vy_data;
bias_data(:,2)=Vg_data;
Vb=0;
bias_data(:,3)=Vb;
Vx=0;
bias_data(:,4)=Vx;

%% now we run optimization. Make a matrix of initial guess
options = optimset('Display','iter','TolFun',1e-14,'MaxFunEvals',8000,'MaxIter',1500);

lb=[1;1;0;1;0;0.1;50;0.2]; % lower bound constraints
ub=[500;500;0.5;2;0.5;10;1000;0.8]; % upper bound constraints

% Optimization routine 
[coeff_op_tran,resnorm,residual,exitflag] = lsqcurvefit(@(coeff,bias_data)mvs_si_1_1_0(input_params,coeff,bias_data),coeff_init,bias_data,log10(Id_data),lb,ub,options); 
end
