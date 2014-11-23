function comp_inv
  Mn = load('snmos');
  Mp = load('spmos');
  disp('generating models for inverters...');
  comp_inv_help(Mn,Mp,2);
end

function comp_inv_help(Mn,Mp,r)
% This function computes simulation model for inverter circuit based on models of sNMOS and sPMOS.  
%  r is the ratio of the PMOS to NMOS, it is 2 by default. 
  if(nargin<3||isempty(r)), r = 2; end

  assert(strcmp(Mn.META.fab,Mp.META.fab)); assert(strcmp(Mn.META.process,Mp.META.process)); 
  assert(Mn.META.gnd==Mp.META.gnd); assert(Mn.META.vdd==Mp.META.vdd);
  assert(Mn.SIZE.wid==Mp.SIZE.wid); assert(Mn.SIZE.len==Mp.SIZE.len);
  assert(all(Mn.GRID.v0==Mp.GRID.v0)); assert(all(Mn.GRID.nv==Mp.GRID.nv)); assert(all(Mn.GRID.dv==Mp.GRID.dv));
  
  GRID = Mn.GRID;  % in is the source of Mn/Mp, out is the drain of Mn/Mp
  SIZE = Mn.SIZE;  % use the SIZE of snmos
  META = Mn.META; 
  META.name = ['inv_',num2str(r)];
  META.nodes = {'i','o'};
  META.desc = ['Inverter circuit  with p/n ratio',num2str(r)]; 
  
  % data
  in = Mn.data; ip = Mp.data;
  data = -(in+ip*r); % output current
  err = Mn.err+Mp.err*r;
  
  % save file
  save(META.name,'data','err','GRID','SIZE','META');
end
