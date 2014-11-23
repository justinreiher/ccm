 %CCM_OPEN: Initialization CCM 
 %  ccm_open;
function ccm_open(debug)
  if(nargin<1||isempty(debug))
    debug = 0;
  end

  disp('Starting CCM ......'); 

  % error and warnings
  if(debug)
    dbstop if error;
  else
    warning off all;
  end
  

  % Add path
  ccm_addpath;

  disp('CCM initialization complete!');
end %function ccm_open


function ccm_addpath
  ccm_home = ccm_info('ccm_home'); 
  ccm_dirs = ccm_info('ccm_dirs'); 
  disp('Add CCM directories into Matlab search path');
  for i=1:length(ccm_dirs)
    dirname = [ccm_home,'/',ccm_dirs{i}];
    addpath(dirname); 
  end
  disp('Add CCM libraries into Matlab search path');
  ccm_libs = ccm_cfg('get','libRoot');
  addpath(genpath(ccm_libs)); % add all subdirs recursively

  addpath(ccm_home); 
end %function ccm_addpath
