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
            
            if(strcmp(deviceType,'pmos'))
                txParams = this.pmosParams;
            end
            
            if(strcmp(deviceType,'nmos'))
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
                if(strcmp(deviceType,'nmos'))
                    Vin(3,indAuto) = 0;
                end
                if(strcmp(deviceType,'pmos'))
                    Vin(3,indAuto) = Vdd;
                end
            end
                
            widMos = reshape(repmat(widVec,1,numPoints)',numPoints*numDevices,1);
            rlenMos = reshape(repmat(rlenVec,1,numPoints)',numPoints*numDevices,1);
            
            txParams.W = widMos';
            txParams.Lgdr = txParams.Lgdr.*rlenMos';
            
            [Ids,capInfo] = mvs_Id(txParams,Vin);
            factor = repmat(I_factor,1,numPoints*numDevices);
            
            
            I = repmat(Ids,numTerminals,1);
            I = I.*factor;
            
            Imos = reshape(I,numTerminals*numDevices,[]);
        end
        
        function Cmos = CmosEval(this,deviceParams,deviceType,Vin,tbOptions)
            
      %get the capacitance model to use
            capModel = tbOptions.capModel;
            
            if(strcmp(deviceType,'pmos'))
                txParams = this.pmosParams;
            end
            
            if(strcmp(deviceType,'nmos'))
                txParams = this.nmosParams;
            end
            
            [numTerminals,numPoints,numDevices] = size(Vin);
            V  = reshape(Vin,numTerminals,[]);
            
            switch(capModel)
                case 'gnd'
                    
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
                    
                case 'full'
                    
                    try
                        gradientinit(1);
                    catch
                        error('IntLAB AD package required for this functionality - get at:  http://www.ti3.tu-harburg.de/intlab/')
                    end
                    
                    devMap = deviceParams.deviceMap;
                    numNodes = tbOptions.numNodes;
                    
                    %%%% Place holder for now
                    widVec = [deviceParams.wid{:}];
                    rlenVec = [deviceParams.rlen{:}];
                    %V = reshape(Vin,numTerminals-1,[]);
                    widMos = reshape(repmat(widVec,1,numPoints)',numPoints*numDevices,1);
                    rlenMos = reshape(repmat(rlenVec,1,numPoints)',numPoints*numDevices,1);
                    devMap  = repmat(devMap,1,numPoints);
                    %%%%
                    
                    jCap = this.junctionCap;
                    jCap.W = widMos';
                    jCap.Lgdr = (rlenMos.*jCap.Lgdr)';
                    
                    Cmos = zeros(numDevices*numTerminals*numPoints,numNodes);
                    
                    [Ipre,~] = mvs_Id(txParams,V);
                    
                    Vgrad = gradientinit(V);
                    [~,capGrad] = mvs_Id_AD(txParams,Vgrad,Ipre);
                    deviceParams.capParams = capGrad;
                    
                    deviceParams.capParams.junctionCap = jCap;
                    
                    Cdev = mvs_c(txParams,Vgrad,deviceParams.capParams);
                    
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
                    
                otherwise
                    error(strcat(capModel,' is an unknown Capacitor model'));
                    
            end
        
        
        end
        
            
        function [residualCurrent,jacobian] = devRes(this,deviceParams,Vin,I,deviceType,autoConnect,numCircuitNodes,tbOptions)
            
            try
                gradientinit(1);
            catch
                error('IntLAB AD package required for this functionality - get at:  http://www.ti3.tu-harburg.de/intlab/')
            end
            
            widVec  = [deviceParams.wid{:}];
            rlenVec = [deviceParams.rlen{:}];
            
            orgSign = I(1,:);
            ind = find(orgSign == -1);
            I = I(2,:);
            
            Vdd = this.modelParams.vdd;
            
            if(strcmp(deviceType,'pmos'))
                txParams = this.pmosParams;
            end
            
            if(strcmp(deviceType,'nmos'))
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
    end
    
end

