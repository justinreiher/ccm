% This function computes simulation model for a transistor whose source is 
% connected to substrate, i.e., ground for NMOS and power for PMOS
% "s" is for substrate
function [Msn, Msp] = comp_smos(Mn,Mp,path)
  if(nargin<3||isempty(path)), path = '.'; end;

  disp('processing nmos ...');
  Msn = comp_smos_help(1,Mn,path); 
  disp('processing pmos ...');
  Msp = comp_smos_help(2,Mp,path);
end
  
function model = comp_smos_help(nORp,M,path)
  %M = ccm_getModel(device);
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
 
  file = [path,'/',META.name,'.mat']; 
  save(file,'data','err','GRID','SIZE','META');
  model = load(file);
end
