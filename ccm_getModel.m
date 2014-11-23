function M = ccm_getModel(name) 
% M = ccm_getModel(name) 
% The function returns model data of specified device.
% The location of the mat file is calculated as 
%   <matRoot>/@matFileFunc(struct(<fab>,<process>,<type>,<name>)).mat
% M has the same contents with the mat file, which must contains: 
%   data: model data for each grid 
%   err:  error term for each grid, could be [] or scalar
%   GRID: grid information
%     v0: initial value for grid 0
%     nv: number of grids
%     dv: unit length of grids
%   SIZE: circuit size
%     len,wid: length and width
%   META: meta-info
%     type: 'simu' or 'quad' 
% 
% For the default mat files, the following information are provided under 'META'
%     lib: library name
%     name: unique name of the device
%     fab:  fab, e.g. TSMC, INTEL, PTM 
%     process: e.g. 180nm, 10nm. 
%     gnd/vdd: voltage for gnd/vdd
%     nodes: nodes names
%     desc: description of the circuit
%
% For 'simu' type, pre-compute the sum of the table for least square method.
%   M.YC: accumulated sum of data for 'simu' type

% Find the mat file to be loaded.
matRoot = ccm_cfg('get','matRoot');
matFileFunc = ccm_cfg('get','matFileFunc');
device.fab = ccm_cfg('get','fab');
device.process = ccm_cfg('get','process'); 
device.type = ccm_cfg('get','type'); 
device.name = name;

% get it directly if loaed before
file = [matRoot,'/',matFileFunc(device)]; 
M = utils_hashtable('get',file);

% load from mat file
if(isempty(M))
  % load the mat file, result is a structure
  disp(['Loading ',file,' ...']);
  M = load(file);
  % pre-computation for 'simu' type
  if(strcmp(lower(device.type),'simu')) 
    YC = M.data;
    for i=1:length(M.GRID.v0) 
      YC = cumsum(YC,i);
    end
    M.YC = YC;
  end
  % save the result 
  utils_hashtable('insert',file,M);
end
