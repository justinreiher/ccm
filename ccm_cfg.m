% function [val,status] = ccm_cfg(op,varargin)
%   This file keeps a global structure to store info shared over CCM pkgs. 
%   It can be use to config CCM or store global data
%   It supports two operations: 
%     GET: 
%       S = ccm_cfg('get');        % get the whole struct
%       value = ccm_cfg('get',key) % get the value
%       {v1,v2,...} = ccm_cfg('get',k1,k2,...) % get multiple values
%     SET:
%       S = ccm_cfg('set',k1,v1,k2,v2,...);  
% 
%   CCM configs
%     matRoot: root directory for the model files.  
%       default: $CCM_HOME/mat
%     matFileFunc: the function for finding mat files given a device 
%       default: <fab>/<process>/<type>/<name>.mat
%     fab: fabrication information
%       default: 'PTM'
%     process: process information
%       default: '180nm'
%     type: model types
%       values: 'simu','quad'
%       default: 'simu'
%     libRoot: root directory for the libraries. 
%       default: $CCM_HOME/libs/coho
%       NOTE: this must be set before 'ccm_open'. 
%     interpMethod: default interpolation methods for current function
%       values: 'coswin','linear','lookup'
%       default: 'coswin'
%     interpJacMethod: default interpolation methods for Jacobian matrix 
%       values: 'ana'(analytical),'num'(numerical)
%       default: 'ana'
%     interpJacNumMethod: default interpolation methods for numerical Jacobian matrix 
%       values: 'linear','coswin','lookup', 
%       default: 'linear'
%       NOTE: The result  of 'lookup' method is bad, the 'coswin' method for 'simu' models is bad.
%     lftSimuMethod: default linearization methods for Simu models
%       values: 'lls','mm','ls','lp'. 
%       default: 'lls'
%       NOTE: 'ls' and 'lp' are very expensive
%     lftQuadMethod: default linearization methods for Quad models
%       values: 'lip','pt'
%       default: 'lip'
%     
function [val,status] = ccm_cfg(op,varargin)
  % NOTE: Because of the Matlab 2013 version bug, I have to use global vars. 
  %       Please don't modify the value by other functions. 
  %       Persistent vars will be re-inited the first time when path changed.
  %persistent CCM_CFG;
  global CCM_CFG;
  if(isempty(CCM_CFG))
    CCM_CFG= ccm_cfg_default;
    disp('init ccm_cfg');
  end
  [val,status,update] = utils_struct(CCM_CFG,op,varargin{:});
  if(status==0&&update)
    ccm_cfg_check(val);
    CCM_CFG = val; % save the update
  end
end % ccm_cfg;

function cfg = ccm_cfg_default
  ccm_home = ccm_info('ccm_home');
  % global parameters
  fab = 'PTM'; process = '180nm'; type = 'simu'; 
  % the function for finding mat files for a model 
  matRoot = [ccm_home,'/mat'];
  matFileFunc = @(d)([upper(d.fab),'/',lower(d.process),'/',lower(d.type),'/',lower(d.name)]); 
  libRoot = [ccm_home,'/libs/coho'];
  % interp methods
  interpMethod = 'coswin';      % coswin, lookup, linear 
  interpJacMethod = 'ana';         % ana(analystical) num(numerical)
  interpJacNumMethod = 'linear';   % linear,coswin (for quad model only), lookup (bad)
  % linfit method
  lftSimuMethod   = 'lls';         % lls mm, (ls, lp too slow)
  lftQuadMethod   = 'lip';         % lip, pt

  % linfit methods
  cfg = struct('fab',fab,'process',process,'type',type, ...
               'matRoot',matRoot,'matFileFunc',matFileFunc, 'libRoot',libRoot, ...
               'interpMethod', interpMethod, 'interpJacMethod', interpJacMethod, ...
               'lftSimuMethod',lftSimuMethod, 'lftQuadMethod',lftQuadMethod ...
              ); 
end % ccm_cfg_default;

function ccm_cfg_check(val) 
  if(~isa(val.matFileFunc,'function_handle'))
    error('You must set matFileFunc as a function handler');
  end
  % 
  if(any(strcmpi(val.lftSimuMethod,{'ls','lp'})))
    warning('The LS or LP method for linearizing SIMU methods is very expensive.');
  end
end % ccm_cfg_check 
