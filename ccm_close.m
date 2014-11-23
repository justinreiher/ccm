% CCM_CLOSE: Release all CCM resources 
%   ccm_close;
function ccm_close
  ccm_rmpath; % do not remove path
  disp('CCM resources have been released!');
end %function ccm_close

function ccm_rmpath
  ccm_home = ccm_info('ccm_home'); 
  ccm_dirs = ccm_info('ccm_dirs'); 
  ccm_libs = ccm_cfg('get','libRoot');
  disp('Remove CCM libraries from Matlab search path');
  rmpath(genpath(ccm_libs)); % add all subdirs recursively
  disp('Remove CCM directories from Matlab search path');
  for i=1:length(ccm_dirs)
    dirname = [ccm_home,'/',ccm_dirs{i}];
    rmpath(dirname); 
  end
  % rmpath(ccm_home);  % leave the ccm_home to call ccm_open
end %function ccm_rmpath
