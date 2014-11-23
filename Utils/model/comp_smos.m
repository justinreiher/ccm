function comp_smos
% This function computes simulation model for a transistor whose source is 
% connected to substrate, i.e., ground for NMOS and power for PMOS
% "s" is for substrate

  disp('processing nmos ...');
  comp_smos_help(1,'nmos');
  disp('processing pmos ...');
  comp_smos_help(2,'pmos');
end
  
function comp_smos_help(nORp,device)
  M = load(device);
  
  GRID = M.GRID; SIZE = M.SIZE; META = M.META;
  
  GRID.v0 = GRID.v0(2:3); GRID.dv= GRID.dv(2:3); GRID.nv = GRID.nv(2:3); 
  META.name = ['s',META.name];
  META.desc = [META.desc,' with the source connected to substrate'];
  META.nodes = META.nodes(2:3);
  
  if(nORp==1)
    sid = round(grid_v2loc(M.GRID,META.gnd,1));
  else
    sid = round(grid_v2loc(M.GRID,META.vdd,1));
  end
  
  data = reshape(M.data(sid,:,:),reshape(GRID.nv,1,[])); 
  if(isempty(M.err))
    err = [];
  elseif(isscalar(M.err))
    err = M.err;
  else 
    err = reshape(M.err(sid,:,:),reshape(GRID.nv,1,[]));
  end
  
  save(META.name,'data','err','GRID','SIZE','META');
end
