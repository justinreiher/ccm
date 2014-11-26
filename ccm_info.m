
% function val = ccm_info(field)
%   This function returns read-only information for CCM, including: 
%     ccm_home:  root path of CAR software 
%     ccm_dirs:  all CCM dirs to be added in Matlab
%     mat_path:  path for model mat data
%     user:      current user
%     version:   CCM version
%     license:   CCM license
%  Ex: 
%     info = ccm_info;  // return the structure
%     has_cplex = ccm_info('has_cplex'); // has the value
function val = ccm_info(field)
  % NOTE: I use global vars because of the Matlab bug. 
	%       (When the code is in linked dir, persistent vars are re-inited 
	%        when firstly changing to a new directory). 
  %       Please don't modify the value by other functions. 
  %persistent  CCM_INFO;
  global CCM_INFO;
  if(isempty(CCM_INFO)) 
    CCM_INFO = ccm_info_init; % evaluate once
  end
  if(nargin<1||isempty(field))
    val = CCM_INFO;
  else
    val = CCM_INFO.(field);
  end; 
end
function  info = ccm_info_init
  % CCM root path
  ccm_home='/ubc/cs/home/c/chaoyan/ccm'; 

  % CCM directories 
  ccm_dirs = {
    'MSpice',
    'MSpice/libs',
    'MSpice/libs/vs',
    'MSpice/libs/coho',
    'Linfit',
    'Linfit/Simu',
    'Linfit/Quad',
    'Linfit/Annulus',
    'Linfit/Utils',
    'Interp',
    'Model',
    'Utils'};

  % current user
  [~,user] = unix('whoami');
  user = user(1:end-1);

  version = 1.0;
  license = 'bsd';
  
  info = struct('version',version, 'license',license, ...
								'ccm_home',ccm_home, 'ccm_dirs',{ccm_dirs}, 'user',user); 
end % ccm_info

