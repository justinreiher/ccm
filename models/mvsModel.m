classdef mvsModel < model
    %mvsModel is the class model that has all the methods to compute
    %currents for the MVS model from: 
    
    properties
        pmosParams = [];
        nmosParams = [];
        junctionCap = [];
        capPerUnit = [];
    end
    
    methods
        function this = mvsModel(modelName,modelConfig)
            this = this@model(modelName,modelConfig);
            paramHandle = modelConfig.modelParams;
            this.pmosParams = paramHandle('p');
            this.nmosParams = paramHandle('n');
            this.capPerUnit = paramHandle('baseCap');
            this.junctionCap = paramHandle('junc');
            this.configureDefaults(modelConfig);
        end
        
        %evaluate PMOS device current
        function [i,cParams] = iPMOS(this,deviceParams,Vin,tbOptions)
            %try getting the parameter of auto connect, if it doesn't exist
            %then the circuit was build with no auto-connect features.
            try
                vddConnectFlag = [deviceParams.vddFlag{:}];
                [Ipmos,cParams] = this.ImosEval('pmos',deviceParams,vddConnectFlag,Vin);
            catch   
                [Ipmos,cParams] = this.ImosEval('pmos',deviceParams,false,Vin);
            end
            i = Ipmos;
        end
            
        %evaluate NMOS device current
        function [i,cParams] = iNMOS(this,deviceParams,Vin,tbOptions)
            %try getting the parameter of auto connect, if it doesn't exist
            %then the circuit was build with no auto-connect features.
            try
                gndConnectFlag = [deviceParams.gndFlag{:}];
                [Inmos,cParams] = this.ImosEval('nmos',deviceParams,gndConnectFlag,Vin);
            catch
                [Inmos,cParams] = this.ImosEval('nmos',deviceParams,false,Vin);
            end
            i = Inmos;
        end
        
       %evaluate PMOS device Capacitance: Model only considers capacitance
       %to ground at the moment.
        function c = cPMOS(this,deviceParams,Vin,tbOptions)  
            Cpmos = this.CmosEval(deviceParams,'pmos',Vin,tbOptions);          
            c = Cpmos;
        end
            
        %evaluate NMOS device Capacitance: Model only consides capacitance
        %to ground at the moment.
        function c = cNMOS(this,deviceParams,Vin,tbOptions)       
            Cnmos = this.CmosEval(deviceParams,'nmos',Vin,tbOptions);
            c = Cnmos;
        end
        
        function [Imos,capInfo] = ImosEval(this,deviceType,deviceParams,autoConnect,Vin)
            widVec  = [deviceParams.wid{:}];
            rlenVec = [deviceParams.rlen{:}];
            
            Vdd = this.modelParams.vdd;
            
            if(deviceType == 'pmos')
                txParams = this.pmosParams;
            end
            
            if(deviceType == 'nmos')
                txParams = this.nmosParams;
            end
            
            [numTerminals,numPoints,numDevices] = size(Vin);
            
            I_factor=[-1;0;1;0]; %this factor nullifies the body terminal of the transistor
            Vin = reshape(Vin,numTerminals,[]);
            
            %check to see if there are any devices that have the autoConnectFlag
            %set
            if(any(autoConnect))
                
                indAuto = find(autoConnect == true);
                %jr: I beleive that this won't work if there are multiple
                %points to evaluate, unsure. Making a note here just so
                %that if strange things happen later on with evaluating
                %multiple phase space points that this is a reminder that
                %the below may be the culprite.
                if(deviceType == 'nmos')
                    Vin(3,indAuto) = 0;
                end
                if(deviceType == 'pmos')
                    Vin(3,indAuto) = Vdd;
                end
            end
                
            widMos = reshape(repmat(widVec,1,numPoints)',numPoints*numDevices,1);
            rlenMos = reshape(repmat(rlenVec,1,numPoints)',numPoints*numDevices,1);
            
            txParams.W = widMos';
            txParams.Lgdr = txParams.Lgdr.*rlenMos';
            
            Ids = mvs_Id(txParams,Vin);
            factor = repmat(I_factor,1,numPoints*numDevices);
            
            
            I = repmat(Ids,numTerminals,1);
            I = I.*factor;
            
            Imos = reshape(I,numTerminals*numDevices,[]);
            capInfo = [];
        end
        
        function Cmos = CmosEval(this,deviceParams,deviceType,Vin,tbOptions)
            
            %get the capacitance model to use
            capModel = tbOptions.capModel;
            
            if(deviceType == 'pmos')
                txParams = this.pmosParams;
            end
            
            if(deviceType == 'nmos')
                txParams = this.nmosParams;
            end
            
            [numTerminals,numPoints,numDevices] = size(Vin);
            C_factor=[1;1;1;1];
            widVec = [deviceParams.wid{:}];
            rlenVec = [deviceParams.rlen{:}];
            %V = reshape(Vin,numTerminals-1,[]);
            widMos = reshape(repmat(widVec,1,numPoints)',numPoints*numDevices,1);
            rlenMos = reshape(repmat(rlenVec,1,numPoints)',numPoints*numDevices,1);
            cap = this.capPerUnit.*widMos.*rlenMos*txParams.Lgdr;
            %exclude body terminal, hacky fix here!!!!!
            cap = repmat(cap',numTerminals,1);
            factor = repmat(C_factor,1,numPoints*numDevices);
            cap = cap.*factor;
            %hack fix to re-introduce the body terminal....
            %cap = [cap;zeros(1,numPoints*numDevices)];
            
            Cmos = reshape(cap,numTerminals*numDevices,[]);
        
        
        end
        
        function configureDefaults(this,modelConfig)
            %right now all this does is sets the default temperature
            % with more model configuration options maybe change this into
            % a switch statement, or a for loop.
            this.pmosParams.Tjun = modelConfig.temp;
            this.nmosParams.Tjun = modelConfig.temp;
            this.junctionCap.Tjun = modelConfig.temp;
        end
    end
    
end

