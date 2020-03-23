classdef ekvModel < model
    %mvsModel is the class model that has all the methods to compute
    %currents for the MVS model from: 
    
    properties
        pmosParams = [];
        nmosParams = [];
        junctionCap = [];
        capPerUnit = [];
        capScale = [];
    end
    
    methods
        function this = ekvModel(modelName,modelConfig)
            this = this@model(modelName,modelConfig);
            paramHandle = modelConfig.modelParams;
            this.pmosParams = paramHandle('p');
            this.nmosParams = paramHandle('n');
            this.capPerUnit = paramHandle('baseCap');
            this.capScale   = paramHandle('capScale');
            this.junctionCap = paramHandle('junc');
            this.configureDefaults(modelConfig);
        end
        
        %evaluate PMOS device current
        function [i,cParam] = iPMOS(this,deviceParams,Vin,tbOptions)
            %try getting the parameter of auto connect, if it doesn't exist
            %then the circuit was build with no auto-connect features.
            try
                vddConnectFlag = [deviceParams.vddFlag{:}];
                [Ipmos,capInfo] = this.ImosEval('pmos',deviceParams,vddConnectFlag,Vin);
            catch   
                [Ipmos,capInfo] = this.ImosEval('pmos',deviceParams,false,Vin);
            end
            i = Ipmos;
            cParam = capInfo;
        end
            
        %evaluate NMOS device current
        function [i,cParam] = iNMOS(this,deviceParams,Vin,tbOptions)
            %try getting the parameter of auto connect, if it doesn't exist
            %then the circuit was build with no auto-connect features.
            try
                gndConnectFlag = [deviceParams.gndFlag{:}];
                [Inmos,capInfo] = this.ImosEval('nmos',deviceParams,gndConnectFlag,Vin);
            catch
                [Inmos,capInfo] = this.ImosEval('nmos',deviceParams,false,Vin);
            end
            i = Inmos;
            cParam = capInfo;
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
        
        function [Imos,capInfo] = ImosEval(this,deviceType,deviceParams,autoConnect,Vin)
            
            %This is to not distinguish between AD version of model
            %evaluation and the non AD version. This model has them as the
            %same. Not the best solution here...
            if(iscell(Vin))
                Vin = Vin{1};
            end
            
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
            if(length(sizeRes) == 2)
                numDevices = 1;
            else
                numDevices   = sizeRes(3);
            end
            
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
            txParams.Lg = txParams.Lg.*rlenMos';
            
            [Ids,i_f,i_r,n] = ekv_Ids(txParams,Vin);
            factor = repmat(I_factor,1,numPoints*numDevices);
            
            
            I = repmat(Ids,numTerminals,1);
            I = I.*factor;
            capInfo.i_f = i_f;
            capInfo.i_r = i_r;
            capInfo.n   = n;
            Imos = reshape(I,numTerminals*numDevices,[]);
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
            
            %broken out like this because hessian/gradient variables do not
            %broadcast the size
            A = size(Vin);
            numTerminals = A(1);
            numPoints    = A(2);
            if(length(A) == 2)
                numDevices = 1;
            else
                numDevices   = A(3);
            end
            V  = reshape(Vin,numTerminals,[]);
            
            switch(capModel)
                case 'gnd'
            
                    C_factor=[1;1;1;1];
                    widVec = [deviceParams.wid{:}];
                    rlenVec = [deviceParams.rlen{:}];
                    %V = reshape(Vin,numTerminals-1,[]);
                    widMos = reshape(repmat(widVec,1,numPoints)',numPoints*numDevices,1);
                    rlenMos = reshape(repmat(rlenVec,1,numPoints)',numPoints*numDevices,1);
                    cap = this.capPerUnit.*widMos.*rlenMos*txParams.Lg;
                    %exclude body terminal, hacky fix here!!!!!
                    cap = repmat(cap',numTerminals,1);
                    factor = repmat(C_factor,1,numPoints*numDevices);
                    cap = cap.*factor;
                    %hack fix to re-introduce the body terminal....
                    %cap = [cap;zeros(1,numPoints*numDevices)];
                    
                    Cmos = reshape(cap,numTerminals*numDevices,[]);
                 
                case 'full'
                    
                    devMap = deviceParams.deviceMap;
                    numNodes = tbOptions.numNodes;


                    
                    %%%% Place holder for now
                    widVec = [deviceParams.wid{:}];
                    rlenVec = [deviceParams.rlen{:}];
                    %V = reshape(Vin,numTerminals-1,[]);
                    widMos = reshape(repmat(widVec,1,numPoints)',numPoints*numDevices,1);
                    rlenMos = reshape(repmat(rlenVec,1,numPoints)',numPoints*numDevices,1);
                    devMap  = repmat(devMap,1,numPoints);
                    cap = (this.capPerUnit.*widMos.*rlenMos*txParams.Lg)';
                    %%%%
                    
                    jCap = this.junctionCap;
                    jCap.W = widMos';
                    jCap.Lgdr = (rlenMos.*jCap.Lgdr)';
                    deviceParams.capParams.junctionCap = jCap;
                    deviceParams.capParams.cap = cap;
                    Cmos = zeros(numDevices*numTerminals*numPoints,numNodes)*V(1);
                    Cdev = ekvModelCap(V,deviceParams.capParams);
                    for i = 1:numDevices*numPoints
                        devCap = Cdev(:,:,i);
                        %Capacitances associated with computing
                        %Vdrain_dot
                        Cmos(1+(i-1)*numTerminals,devMap(1,i)) = devCap(1,2) + devCap(1,4) + devCap(1,3);
                        Cmos(1+(i-1)*numTerminals,devMap(2,i)) = -devCap(2,1);
                        Cmos(1+(i-1)*numTerminals,devMap(3,i)) = -devCap(3,1);
                        Cmos(1+(i-1)*numTerminals,devMap(4,i)) = -devCap(4,1);
                        %Capacitances associated with computing
                        %Vgate_dot
                        Cmos(2+(i-1)*numTerminals,devMap(1,i)) = -devCap(1,2);
                        Cmos(2+(i-1)*numTerminals,devMap(2,i)) = devCap(2,1) + devCap(2,3) + devCap(2,4);
                        Cmos(2+(i-1)*numTerminals,devMap(3,i)) = -devCap(3,2);
                        Cmos(2+(i-1)*numTerminals,devMap(4,i)) = -devCap(4,2);
                        %Capacitances associated with computing
                        %Vsource_dot
                        Cmos(3+(i-1)*numTerminals,devMap(1,i)) = -devCap(1,3);
                        Cmos(3+(i-1)*numTerminals,devMap(2,i)) = -devCap(2,3);
                        Cmos(3+(i-1)*numTerminals,devMap(3,i)) = devCap(3,1) + devCap(3,2) + devCap(3,4);
                        Cmos(3+(i-1)*numTerminals,devMap(4,i)) = -devCap(4,3);
                        %Capacitances associated with computing
                        %Vbody_dot
                        Cmos(4+(i-1)*numTerminals,devMap(1,i)) = -devCap(1,4);
                        Cmos(4+(i-1)*numTerminals,devMap(2,i)) = -devCap(2,4);
                        Cmos(4+(i-1)*numTerminals,devMap(3,i)) = -devCap(3,4);
                        Cmos(4+(i-1)*numTerminals,devMap(4,i)) = devCap(4,1)+devCap(4,2) + devCap(4,3);
                    end
                    Cmos = this.capScale*Cmos;
            end
        
        
        end
        
        function [residualCurrent,dIdx] = devRes(this,deviceParams,Vin,I,deviceType,autoConnect,numCircuitNodes,tbOptions)
            widVec  = [deviceParams.wid{:}];
            rlenVec = [deviceParams.rlen{:}];
            
            orgSign = I(1,:);
            ind = find(orgSign == -1);
            I = I(2,:);
            I = orgSign.*I; % no need to strip away sign in this model
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
            
            I_AD = ekv_Ids(txParams,Vin);
            dIdx_Dev = I_AD(:).dx;
            dIdx_Dev = reshape(dIdx_Dev,numDevices,[]);
            dIdV = dIdx_Dev(:,1:numCircuitNodes);
            %         factor = repmat(I_factor,1,numPoints*numDevices);
            %        Jfactor = repmat(J_factor,1,numPoints*numDevices);
            
            
            dIdx = zeros(numDevices*numTerminals,numCircuitNodes);
            for i = 1:numDevices
                dIdx(1+(i-1)*numTerminals:i*numTerminals,:) = I_factor.*dIdV(i,:);
            end
            
            I_AD = repmat(I_AD,numTerminals,1);
            I_AD = I_AD.x;
            I_AD = I_factor.*I_AD;
            
            residualCurrent = I_AD-I_factor.*I;
           
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
    end
    
end

