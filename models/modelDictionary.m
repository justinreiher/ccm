classdef modelDictionary
    %jr: Version 1.0 of the model dictionary contains:
    % model inter: an interpolated model currently fit for PTM 180nm process
    % model MVS  : a short channel model currently fit for the PTM 45nm
    % process
    % model Long Channel : a simple long channel model
    %to add additional models simply add the name of the model to the
    %case statement and implement the required classes
    properties
        version = '1.0';
    end
    
    methods
        function this = modelDictionary()
            %sets no properties. Creates the class provides a method to get
            %specific models
        end
        
        function model = getModel(this,myModel,modelProcess,tbOptions)
            %function that takes in a model name and dispatches the correct
            %model
            
            %if it is referred to as a keyword then load the appropriate
            %configuration. Can be a struct with all the necessary information
            %already provided in which case it is assumed to be correct,
            %otherwise the modelConfig key is used to get the specific
            %configuration.
            if(~isstruct(simOptions))
                modelProcess = this.getModelConfig(myModel,modelProcess,tbOptions);
            end
            
            switch(myModel)
                case 'interp'
                    model = interpModel(myModel,modelProcess);
                case 'MVS'
                    model = mvsModel(myModel,modelProcess);
                case 'MVS AD'
                    model = mvsADModel(myModel,modelProcess);
                case 'EKV'
                    model = ekvModel(myModel,modelProcess);
                case 'EKV AD'
                    model = ekvModel(myModel,modelProcess); %no difference between AD version and non-AD version
                otherwise
                    error(strcat(myModel,' is not a supported model yet'))
            end
        end
        
        % selects the right configuration function to configure the model
        % with the applied testbench options
        
        function modelConfig = getModelConfig(this,model,config,tbOptions)
            
            switch(model)
                case 'interp'
                    modelConfig = ccm_cgf_interp(config,tbOptions);
                case 'MVS'
                    modelConfig = ccm_cgf_mvs(config,tbOptions);
                case 'MVS AD'
                    modelConfig = ccm_cgf_mvsAD(config,tbOptions);
                case 'EKV'
                    modelConfig = ccm_cgf_ekv(config,tbOptions);
                case 'EKV AD'
                    modelConfig = ccm_cgf_ekv(config,tbOptions);
                otherwise
                    error(strcat(config,' is not a valid configuration key'));
            end
        end
    end
end

