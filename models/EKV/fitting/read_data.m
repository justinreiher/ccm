% Written by Yan Peng, UBC, 2016/10/19
% Modified by Justin Reiher, UBC, 2017/05/30
% Re-written by Justin Reiher, UBC, 2019/03/27 - avoid pmos and nmos
% confusion

function [ids_sweepVd_N, ids_sweepVg_N, ids_sweepVs_N,...
    ids_sweepVd_P, ids_sweepVg_P, ids_sweepVs_P] = read_data(nmosSweepVd, nmosSweepVg,nmosSweepVs,...
                                                                      pmosSweepVd, pmosSweepVg,pmosSweepVs, sweeps)
%% Read the nmos data first. Vs = 0 in this data

Vmax = sweeps(1);
%fine grain sweeping
Vx_start  = sweeps(2);
Vx_step_size   = sweeps(3);
Vx_end       = sweeps(4);

%course grain sweeping
Vy_start     = sweeps(5);
Vy_step_size = sweeps(6);
Vy_end       = sweeps(7);

%%% i.e Sweep Vg in step sizes of 0.1 and Vd in step sizes of 0.01

Vx = Vx_start:Vx_step_size:Vx_end;
Vy = Vy_start:Vy_step_size:Vy_end;


%% read in nmos data
%format of Vd sweep file: Vd|Vd|Vg|Ids
idvd_n_raw = dlmread(nmosSweepVd);

%format of Vg sweep file: Vg|Vd|Vg|Ids
idvg_n_raw = dlmread(nmosSweepVg);

%format of Vs sweep file: Vs|Vg|Vs|Ids
idvs_n_raw = dlmread(nmosSweepVs);

% To avoid confusion save Vd|Vg|Vs|Vb|Ids

ids_sweepVd_N = zeros(length(Vx)*length(Vy),5);
ids_sweepVg_N = zeros(length(Vx)*length(Vy),5);
ids_sweepVs_N = zeros(length(Vx)*length(Vy),5);


%(body is tied to zero in the n-device)

for i = 1:length(Vy)
    %popoluate Vd
    ids_sweepVd_N(1+(i-1)*length(Vx):i*length(Vx),1) = idvd_n_raw(1+(i-1)*length(Vx):i*length(Vx),2);
    ids_sweepVg_N(1+(i-1)*length(Vx):i*length(Vx),1) = idvg_n_raw(1+(i-1)*length(Vx):i*length(Vx),2);
    ids_sweepVs_N(1+(i-1)*length(Vx):i*length(Vx),1) = zeros(size(Vx));
    
    %populate Vg
    ids_sweepVd_N(1+(i-1)*length(Vx):i*length(Vx),2) = idvd_n_raw(1+(i-1)*length(Vx):i*length(Vx),3);
    ids_sweepVg_N(1+(i-1)*length(Vx):i*length(Vx),2) = idvg_n_raw(1+(i-1)*length(Vx):i*length(Vx),3);
    ids_sweepVs_N(1+(i-1)*length(Vx):i*length(Vx),2) = idvs_n_raw(1+(i-1)*length(Vx):i*length(Vx),2);
    
    %populate Vs
    ids_sweepVd_N(1+(i-1)*length(Vx):i*length(Vx),3) = zeros(size(Vx));
    ids_sweepVg_N(1+(i-1)*length(Vx):i*length(Vx),3) = zeros(size(Vx));
    ids_sweepVs_N(1+(i-1)*length(Vx):i*length(Vx),3) = idvs_n_raw(1+(i-1)*length(Vx):i*length(Vx),3);
    
    %populate Vb
    ids_sweepVd_N(1+(i-1)*length(Vx):i*length(Vx),4) = zeros(size(Vx));
    ids_sweepVg_N(1+(i-1)*length(Vx):i*length(Vx),4) = zeros(size(Vx));
    ids_sweepVs_N(1+(i-1)*length(Vx):i*length(Vx),4) = zeros(size(Vx));
    
    %populate Ids
    ids_sweepVd_N(1+(i-1)*length(Vx):i*length(Vx),5) = idvd_n_raw(1+(i-1)*length(Vx):i*length(Vx),4);
    ids_sweepVg_N(1+(i-1)*length(Vx):i*length(Vx),5) = idvg_n_raw(1+(i-1)*length(Vx):i*length(Vx),4);
    ids_sweepVs_N(1+(i-1)*length(Vx):i*length(Vx),5) = idvs_n_raw(1+(i-1)*length(Vx):i*length(Vx),4);
    

end

%% read in ppmos data
%format of Vd sweep file: Vd|Vd|Vg|Vs|Ids
idvd_p_raw = dlmread(pmosSweepVd);

%format of Vg sweep file: Vg|Vd|Vg|Vs|Ids
idvg_p_raw = dlmread(pmosSweepVg);

%format of Vs sweep file: Vs|Vd|Vg|Vs|Ids
idvs_p_raw = dlmread(pmosSweepVs);

% To avoid confusion save Vd|Vg|Vs|Vb|Ids

ids_sweepVd_P = zeros(length(Vx)*length(Vy),5);
ids_sweepVg_P = zeros(length(Vx)*length(Vy),5);
ids_sweepVs_P = zeros(length(Vx)*length(Vy),5);


%(body is tied to Vmax in the p-device)

for i = 1:length(Vy)
    %popoluate Vd
    ids_sweepVd_P(1+(i-1)*length(Vx):i*length(Vx),1) = idvd_p_raw(1+(i-1)*length(Vx):i*length(Vx),2);
    ids_sweepVg_P(1+(i-1)*length(Vx):i*length(Vx),1) = idvg_p_raw(1+(i-1)*length(Vx):i*length(Vx),2);
    ids_sweepVs_P(1+(i-1)*length(Vx):i*length(Vx),1) = Vmax*ones(size(Vx));
    
    %populate Vg
    ids_sweepVd_P(1+(i-1)*length(Vx):i*length(Vx),2) = idvd_p_raw(1+(i-1)*length(Vx):i*length(Vx),3);
    ids_sweepVg_P(1+(i-1)*length(Vx):i*length(Vx),2) = idvg_p_raw(1+(i-1)*length(Vx):i*length(Vx),3);
    ids_sweepVs_P(1+(i-1)*length(Vx):i*length(Vx),2) = idvs_p_raw(1+(i-1)*length(Vx):i*length(Vx),3);
    
    %populate Vs
    ids_sweepVd_P(1+(i-1)*length(Vx):i*length(Vx),3) = Vmax*ones(size(Vx));
    ids_sweepVg_P(1+(i-1)*length(Vx):i*length(Vx),3) = Vmax*ones(size(Vx));
    ids_sweepVs_P(1+(i-1)*length(Vx):i*length(Vx),3) = idvs_p_raw(1+(i-1)*length(Vx):i*length(Vx),4);
    
    %populate Vb
    ids_sweepVd_P(1+(i-1)*length(Vx):i*length(Vx),4) = Vmax*ones(size(Vx));
    ids_sweepVg_P(1+(i-1)*length(Vx):i*length(Vx),4) = Vmax*ones(size(Vx));
    ids_sweepVs_P(1+(i-1)*length(Vx):i*length(Vx),4) = Vmax*ones(size(Vx));
    
    %populate Ids
    ids_sweepVd_P(1+(i-1)*length(Vx):i*length(Vx),5) = idvd_p_raw(1+(i-1)*length(Vx):i*length(Vx),5);
    ids_sweepVg_P(1+(i-1)*length(Vx):i*length(Vx),5) = idvg_p_raw(1+(i-1)*length(Vx):i*length(Vx),5);
    ids_sweepVs_P(1+(i-1)*length(Vx):i*length(Vx),5) = idvs_p_raw(1+(i-1)*length(Vx):i*length(Vx),5);
    

end

end
