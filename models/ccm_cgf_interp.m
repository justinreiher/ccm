function cgf = ccm_cgf_interp(key,simOptions)
%ccm_cgf_interp returns a configuration struct for the interpolated model
%used in the simulator.
%To add additional models, simply add it to the switch statement and place
%the appropriate mat file in the appropriate location and it should be able
%to find the relevant information.
warning('simOptions are not supported here')
switch(key)
    case 'PTM 180nm'
        cgf = ccm_cgf('PTM','180nm');
    otherwise
        error(strcat(key,' unknown configuration key'));
end


end

function config = ccm_cgf(fab,process)
ccm_home = ccm_info('ccm_home');
  % global parameters
  type = 'simu'; 
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
  config = struct('fab',fab,'process',process,'type',type, ...
               'matRoot',matRoot,'matFileFunc',matFileFunc, 'libRoot',libRoot, ...
               'interpMethod', interpMethod, 'interpJacMethod', interpJacMethod, ...
               'lftSimuMethod',lftSimuMethod, 'lftQuadMethod',lftQuadMethod ...
              ); 
end 