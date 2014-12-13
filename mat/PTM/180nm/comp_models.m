% This script is used to compute models for snmos/spmos/inv based on nmos/pmos models
% Also compute quadratic polynomial models based on simulation models

cd simu
disp('Generate models for inverters and nmos/pmos with source connected to ground/vdd');
Mn = load('nmos.mat');
Mp = load('pmos.mat');
comp_smos(Mn,Mp,'.');
comp_inv(Mn,Mp,2,'.');

% comp_quadModels;
disp('Generate quadratic models');
grid.v0 = -0.2; grid.dv = 0.1; grid.nv = 21; % [-0.1,1.9]
models = {'inv_2','snmos','spmos','nmos','pmos'};
for i=1:length(models)
  fprintf('Working on %s\n',models{i});
  M = load(models{i});
  comp_quadModel(M,grid,'../quad');
end

