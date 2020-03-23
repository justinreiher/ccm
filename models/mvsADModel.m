classdef mvsADModel < model
    %mvsModel is the class model that has all the methods to compute
    %currents for the MVS model from:
    
    properties
        pmosParams = [];
        nmosParams = [];
        junctionCap = [];
        capPerUnit = [];
    end
    
    methods
        function this = mvsADModel(modelName,modelConfig)
            this = this@model(modelName,modelConfig);
            paramHandle = modelConfig.modelParams;
            this.pmosParams = paramHandle('p');
            this.nmosParams = paramHandle('n');
            this.junctionCap = paramHandle('junc');
            this.capPerUnit = paramHandle('baseCap');
            this.configureDefaults(modelConfig);
        end
        
        %evaluate PMOS device current
        function i = iPMOS(this,deviceParams,Vin,tbOptions)
            %try getting the parameter of auto connect, if it doesn't exist
            %then the circuit was build with no auto-connect features.
            try
                vddConnectFlag = [deviceParams.vddFlag{:}];
                Ipmos = this.ImosEval('pmos',deviceParams,vddConnectFlag,Vin);
            catch
                Ipmos = this.ImosEval('pmos',deviceParams,false,Vin);
            end
            i = Ipmos;
        end
        
        %evaluate NMOS device current
        function i = iNMOS(this,deviceParams,Vin,tbOptions)
            %try getting the parameter of auto connect, if it doesn't exist
            %then the circuit was build with no auto-connect features.
            try
                gndConnectFlag = [deviceParams.gndFlag{:}];
                Inmos = this.ImosEval('nmos',deviceParams,gndConnectFlag,Vin);
            catch
                Inmos = this.ImosEval('nmos',deviceParams,false,Vin);
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
        
        %evaluate residual current in NMOS device
        function [residualCurrent,jacobian] = iresNMOS(this,deviceParams,V,I,numCircuitNodes,tbOptions)
            try
                gndConnectFlag = [deviceParams.gndFlag{:}];
                [ires,J] = this.devRes(deviceParams,V,I,'nmos',gndConnectFlag,numCircuitNodes,tbOptions);
            catch
                [ires,J] = this.devRes(deviceParams,V,I,'nmos',false,numCircuitNodes,tbOptions);
            end
            residualCurrent = ires;
            jacobian = J;
        end
        
        %ievaluate residual current in PMOS device
        function [residualCurrent,jacobian] = iresPMOS(this,deviceParams,V,I,numCircuitNodes,tbOptions)
            try
                vddConnectFlag = [deviceParams.vddFlag{:}];
                [ires,J] = this.devRes(deviceParams,V,I,'pmos',vddConnectFlag,numCircuitNodes,tbOptions);
            catch
                [ires,J] = this.devRes(deviceParams,V,I,'pmos',false,numCircuitNodes,tbOptions);
            end
            residualCurrent = ires;
            jacobian = J;
        end
        
        function Imos = ImosEval(this,deviceType,deviceParams,autoConnect,inputs)
            
            Vin = inputs{1};
            IpreComputed = inputs{2};
            
            widVec  = [deviceParams.wid{:}];
            rlenVec = [deviceParams.rlen{:}];
            
            Vdd = this.modelParams.vdd;
            
            if(deviceType == 'pmos')
                txParams = this.pmosParams;
            end
            
            if(deviceType == 'nmos')
                txParams = this.nmosParams;
            end
            
            %the below is broken out like this because size() of a gradient
            %variable does not allow broadcasting results, i.e.
            %[numTerm,numPoint,numDev] = size(Vin) for whatever reason.
            sizeRes = size(Vin);
            numTerminals = sizeRes(1);
            numPoints    = sizeRes(2);
            numDevices   = sizeRes(3);
            
            
            I_factor=[-1;0;1;0]; %this factor nullifies the body terminal of the transistor
            Vin = reshape(Vin,numTerminals,[]);
            
            %check to see if there are any devices that have the autoConnectFlag
            %set
            if(any(autoConnect))
                
                VinTemp = Vin.x;
                indAuto = find(autoConnect == true);
                
                if(deviceType == 'nmos')
                    VinTemp(3,indAuto) = 0;
                end
                if(deviceType == 'pmos')
                    VinTemp(3,indAuto) = Vdd;
                end
                
                %if the input voltage was setup as a gradient AD variable,
                %re-initialize Vin with the autoconnected voltages:
                if(isa(Vin,'gradient'))
                    Vin = gradientinit(VinTemp);
                end
                
                %if the input voltage was setup as a hessian AD variable,
                %re-initialize Vin with the autoconnected voltages:
                if(isa(Vin,'hessian'))
                    Vin = hessianinit(VinTemp);
                end
            end
            
            %            if(isa(widMos,'gradient')||isa(widMos,'hessian'))
            %                [Vin,widMos] = this.augmentADVariables(Vin,widMos);
            %            end
            
            widMos = reshape(repmat(widVec,1,numPoints)',numPoints*numDevices,1);
            rlenMos = reshape(repmat(rlenVec,1,numPoints)',numPoints*numDevices,1);
            
            txParams.W = widMos';
            txParams.Lgdr = txParams.Lgdr.*rlenMos';
            
            Ids = mvs_Id_AD(txParams,Vin,IpreComputed);
            factor = repmat(I_factor,1,numPoints*numDevices);
            
            
            I = repmat(Ids,numTerminals,1);
            I = I.*factor;
            
            Imos = reshape(I,[numTerminals,numPoints,numDevices]);
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
            if(capModel == 'gnd')
                
                
                
                [numTerminals,numPoints,numDevices] = size(Vin.x);
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
                
                Cmos = reshape(cap,[numTerminals,numPoints,numDevices]);
            else
                error(strcat(capModel,' is an unknown Capacitor model'));
                
            end
            
            
        end
        
        function [residualCurrent,jacobian] = devRes(this,deviceParams,Vin,I,deviceType,autoConnect,numCircuitNodes,tbOptions)
            widVec  = [deviceParams.wid{:}];
            rlenVec = [deviceParams.rlen{:}];
            
            orgSign = I(1,:);
            ind = find(orgSign == -1);
            I = I(2,:);
            
            Vdd = this.modelParams.vdd;
            
            if(deviceType == 'pmos')
                txParams = this.pmosParams;
            end
            
            if(deviceType == 'nmos')
                txParams = this.nmosParams;
            end
            
            %the below is broken out like this because size() of a gradient
            %variable does not allow broadcasting results, i.e.
            %[numTerm,numPoint,numDev] = size(Vin) for whatever reason.
            sizeRes = size(Vin);
            numTerminals = sizeRes(1);
            numPoints    = sizeRes(2);
            numDevices   = sizeRes(3);
            
            
            I_factor=[-1;0;1;0]; %this factor nullifies the body terminal of the transistor
            Vin = reshape(Vin,numTerminals,[]);
            
            %check to see if there are any devices that have the autoConnectFlag
            %set
            if(any(autoConnect))
                
                VinTemp = Vin.x;
                indAuto = find(autoConnect == true);
                
                if(deviceType == 'nmos')
                    VinTemp(3,indAuto) = 0;
                end
                if(deviceType == 'pmos')
                    VinTemp(3,indAuto) = Vdd;
                end
                
                %if the input voltage was setup as a gradient AD variable,
                %re-initialize Vin with the autoconnected voltages:
                if(isa(Vin,'gradient'))
                    Vin = gradientinit(VinTemp);
                end
                
                %if the input voltage was setup as a hessian AD variable,
                %re-initialize Vin with the autoconnected voltages:
                if(isa(Vin,'hessian'))
                    Vin = hessianinit(VinTemp);
                end
            end
            
            widMos = reshape(repmat(widVec,1,numPoints)',numPoints*numDevices,1);
            rlenMos = reshape(repmat(rlenVec,1,numPoints)',numPoints*numDevices,1);
            
            txParams.W = (widMos.*rlenMos)';
            %txParams.Lgdr = txParams.Lgdr*rlenMos';
            VI_AD = gradientinit([Vin.x;I]);
            Ires = computeResidualIdsAD(txParams,VI_AD);
            Jdev = Ires(:).dx;
            Jdev = reshape(Jdev,numDevices,[]);
            Jx = Jdev(:,1:numCircuitNodes);
            Jy = Jdev(:,numCircuitNodes+1:end);
            J = -Jy\Jx;
            %         factor = repmat(I_factor,1,numPoints*numDevices);
            %        Jfactor = repmat(J_factor,1,numPoints*numDevices);
            
            IresSign = sign(Ires.x);
            
            jacobian = zeros(numDevices*numTerminals,numCircuitNodes);
            for i = 1:numDevices
                jacobian(1+(i-1)*numTerminals:i*numTerminals,:) = orgSign(i)*I_factor.*J(i,:);
            end
            
            Ires = repmat(Ires,numTerminals,1);
            Ires = Ires.x;
            Ires = I_factor.*Ires;
            
            residualCurrent = Ires;
           
            %           jacobian = J;
        end
        
        function configureDefaults(this,modelConfig)
            %right now all this does is sets the default temperature
            % with more model configuration options maybe change this into
            % a switch statement, or a for loop.
            this.pmosParams.Tjun = modelConfig.temp;
            this.nmosParams.Tjun = modelConfig.temp;
            this.junctionCap.Tjun = modelConfig.temp;
        end
        
        function [Vin,mosWid] = augmentADVariables(Vin,mosWid)
            tempVin = Vin.x;
            [~,lenVin] = size(tempVin);
            tempMosWid = mosWid.x;
            
            if(lenVin ~= length(tempMosWid))
                tempMosWid = ones(1,lenVin)*tempMosWid;
            end
            
            if(isa(Vin,'hessian')||isa(mosWid,'hessian'))
                augAD = hessianinit([tempVin;tempMosWid]);
            else
                augAD = gradientinit([tempVin;tempMosWid]);
            end
            
            Vin = augAD(1:lenVin,:);
            mosWid = augAD(lenVin+1:end,:);
        end
        
    end
end

