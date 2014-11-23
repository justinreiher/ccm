function spice_simu_ptm
  devices = {'nmos','pmos'}; 
  %devices = {'pmos2','nmos2'};
  for i=1:length(devices)
    fprintf('Working on %s\n',devices{i});
    spice_help(devices{i});
  end
end

function spice_help(device,ids)
  if(nargin<2||isempty(ids))
    ids = true;
  end

  lib = 'coho';
  type = 'simu';
  wid = 10e-6; len = 180e-9; 
  ss = -0.3; ee = 2.25; dv = 0.01; v = ss:dv:ee;
  v0 = ones(3,1)*(ss-dv);
  dv = ones(3,1)*dv; 
  nv = repmat(length(v),3,1);
  gnd=0; vdd=1.8; process=180;
  
  switch(device)
    case 'nmos'
      nodes = {'s','g','d'};
      info = 'PTM 180nm NMOS';
    case 'pmos'
      nodes = {'s','g','d'};
      info = 'PTM 180nm PMOS';
    case 'nmos2'
      nodes = {'g1','g2','d'};
      info = 'PTM 180ns two serial NMOS';
      wid = [wid;wid];
    case 'pmos2'
      nodes = {'g1','g2','d'};
      info = 'PTM 180ns two serial PMOS';
      wid = [wid;wid];
    otherwise
      error('do not support');
  end
  grids = {v;v;v};
  opt.wid = wid;  opt.device = device; opt.ids = ids;
  tmpName = device;
  data = spice_device_ptm(grids,opt,tmpName);
  
  file = [type,'_',device];
  if(~ids)
    file = [file,'_isd'];
  end
  wid = wid(1); % make sure it is a scalar
  save(file,'lib,','type','device','process','gnd','vdd','nodes','info',...
       'wid','len', 'v0','dv','nv','data'); 
end
