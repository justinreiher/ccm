% I found large mismatch between simu and quad models. 
% It is only for nmos, pmos, not for inverter. Debug it.
% check TSMC models to compare
% TSMC has similar problems. 
% UPDATE: This is not a bug. In fact, there is always difference between simu and quad models. 
% For some region, the error is visiablly "big" as current are very tiny. So not a problem.

  name = 'nmos'; 

  ccm_cfg('set','fab','ptm');
  ccm_cfg('set','type','simu');
  m1 = ccm_getModel(name);
  ccm_cfg('set','type','quad');
  m2 = ccm_getModel(name);

  ccm_cfg('set','fab','tsmc');
  ccm_cfg('set','type','simu');
  m3 = ccm_getModel(name);
  ccm_cfg('set','type','quad');
  m4 = ccm_getModel(name);

  grids = m1.META.gnd:0.01:m1.META.vdd;
  V = repmat([0.0;0.5;1.5],1,length(grids));
  V(1,:) = grids;

  ids1 = interp_data(m1,V,[],'coswin');
  ids2 = interp_data(m2,V,[],'coswin');
  ids3 = interp_data(m3,V,[],'coswin');
  ids4 = interp_data(m4,V,[],'coswin');

  figure; hold on;
  plot(grids,ids1,'r');
  plot(grids,ids2,'g');
  plot(grids,ids1,'b');
  plot(grids,ids2,'k');

