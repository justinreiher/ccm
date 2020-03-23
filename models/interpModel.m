classdef interpModel < model
    %interpModel is the class model that has all the methods to compute
    %currents for the interpolated model used in Chao's original work
    
    properties(Constant)
        capPerUnit = 2e-9; %no idea where this number came from...
    end
    
    methods
        function this = interpModel(modelName,modelParameters)
            this = this@model(modelName,modelParameters);
        end
        
        %evaluate PMOS device current
        function i = iPMOS(this,deviceParams,Vin,simOptions)  
            try
                vddConnectFlag = [deviceParams.vddFlag{:}];
                Ipmos = this.ImosEval('pmos',deviceParams,vddConnectFlag,Vin);
            catch
                vddConnectFlag = false;
                Ipmos = this.ImosEval('pmos',deviceParams,vddConnectFlag,Vin);
            end
            i = Ipmos;
        end
            
        %evaluate NMOS device current
        function i = iNMOS(this,deviceParams,Vin,simOptions)
            try
                gndConnectFlag = [deviceParams.gndFlag{:}];
                Inmos = this.ImosEval('nmos',deviceParams,gndConnectFlag,Vin);
            catch
                gndConnectFlag = false;
                Inmos = this.ImosEval('nmos',deviceParams,gndConnectFlag,Vin);
            end
            i = Inmos;
        end
        
       %evaluate PMOS device Capacitance: Model only considers capacitance
       %to ground at the moment.
        function c = cPMOS(this,deviceParams,Vin,simOptions)  
            Cpmos = this.CmosEval(deviceParams,Vin,simOptions);          
            c = Cpmos;
        end
            
        %evaluate NMOS device Capacitance: Model only consides capacitance
        %to ground at the moment.
        function c = cNMOS(this,deviceParams,Vin,simOptions)       
            Cnmos = this.CmosEval(deviceParams,Vin,simOptions);
            c = Cnmos;
        end
        
        function Imos = ImosEval(this,deviceType,deviceParams,autoConnect,Vin)
            widVec  = [deviceParams.wid{:}];
            rlenVec = [deviceParams.rlen{:}];
            
            
            [numTerminals,numPoints,numDevices] = size(Vin);
            
            I_factor=[1;0;-1;0]; %this factor nullifies the body terminal of the transistor
            Vin = reshape(Vin,numTerminals,[]);
            
            %check to see if there are any devices that have the autoConnectFlag
            %set
            I = zeros(numDevices*numPoints);
            if(any(autoConnect))
                indSmos = find(autoConnect == true);
                indMos  = find(autoConnect == false);
                I_factorAutoCon = [-1;0;0;0]; %this factor also nullifies the sources 
                %reshaping of spmos devices width and relative length
                widSmos = widVec(indSmos);
                rlenSmos = rlenVec(indSmos);
                widSmos = reshape(repmat(widSmos,1,numPoints)',numPoints*numDevices,1);
                rlenSmos = reshape(repmat(rlenSmos,1,numPoints)',numPoints*numDevices,1);
                Ismos = interp_device(strcat('s',deviceType),Vin(1:2,:),widSmos,rlenSmos);
                %reshaping of pmos devices width and relative length, only
                %do this if there are devices that are not autoconnected.
                if(~isempty(indMos))
                    widMos = widVec(indMos);
                    rlenMos = rlenVec(indMos);
                    widMos = reshape(repmat(widMos,1,numPoints)',numPoints*numDevices,1);
                    rlenMos = reshape(repmat(rlenMos,1,numPoints)',numPoints*numDevices,1);
                    Imos = interp_device(deviceType,Vin,widMos,rlenMos);
                    %place the current in the appropriate place
                    I(indMos) = Imos;
                end
                %compute the currents         
                I(indSmos) = Ismos;
                %find the columns where the autoConnect factor needs to be
                %placed.
                autoConnect = reshape(repmat(autoConnect,1,numPoints)',numPoints*numDevices,1);
                indFac = autoConnect == true;
                factor = repmat(I_factor,1,numPoints*numDevices);
                factor(:,indFac) = I_factorAutoCon;
            else
                
                widMos = reshape(repmat(widVec,1,numPoints)',numPoints*numDevices,1);
                rlenMos = reshape(repmat(rlenVec,1,numPoints)',numPoints*numDevices,1);           
                I = interp_device(deviceType,Vin(1:3,:),widMos,rlenMos);
                factor = repmat(I_factor,1,numPoints*numDevices);
            end
            
            
            I = repmat(I',numTerminals,1);
            I = I.*factor;
            
            Imos = reshape(I,numTerminals*numDevices,[]);
        end
        
        function Cmos = CmosEval(this,deviceParams,Vin,simOptions)
            [numTerminals,numPoints,numDevices] = size(Vin);
            C_factor=[1;1;1;1];
            widVec = [deviceParams.wid{:}];
            rlenVec = [deviceParams.rlen{:}];
            %V = reshape(Vin,numTerminals-1,[]);
            widMos = reshape(repmat(widVec,1,numPoints)',numPoints*numDevices,1);
            rlenMos = reshape(repmat(rlenVec,1,numPoints)',numPoints*numDevices,1);
            
            cap = this.capPerUnit.*widMos.*rlenMos;
            %exclude body terminal, hacky fix here!!!!!
            cap = repmat(cap',numTerminals,1);
            factor = repmat(C_factor,1,numPoints*numDevices);
            cap = cap.*factor;
            %hack fix to re-introduce the body terminal....
            %cap = [cap;zeros(1,numPoints*numDevices)];
            
            Cmos = reshape(cap,numTerminals*numDevices,[]);
        
        
        end
    end
    
end

