function cgf = ccm_cgf_mvsAD(key,tbOptions)
%ccm_cgf_mvs returns a configuration struct for the mvs model
%used in the simulator with any modifications needed from the simOptions.
%To add additional models, simply add it to the switch statement and place
%the appropriate model_params_<process> function in the appropriate location
%and it should be able to find the relevant information.

switch(key)
    case 'PTM 45nmHP'
        paramHandle = @(type) mvs_model_params_45nmHP(type);
        tbOptions.modelParams = paramHandle;
        cgf = ccm_cgf('PTM','45nmHP',tbOptions);
    otherwise
        error(strcat(key,' unknown configuration key'));
end


end

function config = ccm_cgf(fab,process,tbOptions)
  % global parameters
  config = struct('fab',fab,'process',process); 
  options = fieldnames(tbOptions);
  numFields = length(options);
  
  for i = 1:numFields
      key = options{i};
      config.(key) = tbOptions.(key);
  end

  
end 