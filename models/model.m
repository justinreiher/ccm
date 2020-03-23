classdef model < matlab.mixin.Copyable %handle


    
    properties
        modelName = '';
        modelParams = [];
    end
    
    methods
        function this = model(modelName,modelParams)
            %this dispatches models based on avaialbles models
            this.modelName = modelName;
            this.modelParams = modelParams;
        end   
           
        function [i,cParam] = I(this,deviceType,deviceParams,Vin,simOptions)
            %select the device
            switch(deviceType)
                case 'pmos'
                    [i,cParam] = this.iPMOS(deviceParams,Vin,simOptions);
                case 'nmos'
                    [i,cParam] = this.iNMOS(deviceParams,Vin,simOptions);
                case 'vsrc'
                    i = zeros(size(Vin)); % voltage sources for now produce no current;
                    cParam = []; %by default the voltage source doesn't change the capacitance
                case 'inductor'
                    [i,cParam] = this.iL(deviceParams,Vin,simOptions);
                case 'capacitor'
                    [i,cParam] = this.iC(deviceParams,Vin,simOptions);
                case 'resistor'
                    [i,cParam] = this.iR(deviceParams,Vin,simOptions);
                otherwise
                    error(strcat(deviceType,' is not a supported device type'))
            end
        end
        
        function [ires,J] = IdevRes(this,deviceType,deviceParams,V,I,numCircuitNodes,simOptions)
            %select the device type
            switch(deviceType)
                case 'pmos'
                    [ires,J] = this.iresPMOS(deviceParams,V,I,numCircuitNodes,simOptions);
                case 'nmos'
                    [ires,J] = this.iresNMOS(deviceParams,V,I,numCircuitNodes,simOptions);
                case 'vsrc'
                    [ires,J] = zeros(size(VI)); % voltage sources for now produce no current;
                case 'inductor'
                    [ires,J] = this.iresL(deviceParams,V,I,numCircuitNodes,simOptions);
                case 'capacitor'
                    [ires,J] = this.iresC(deviceParams,V,I,numCircuitNodes,simOptions);
                case 'resistor'
                    [ires,J] = this.iresR(deviceParams,V,I,numCircuitNodes,simOptions);
                otherwise
                    error(strcat(deviceType,' is not a supported device type'))
            end
        end
        
        function c = C(this,deviceType,deviceParams,Vin,simOptions)
            %select the device
            switch(deviceType)
                case 'pmos'
                    c = this.cPMOS(deviceParams,Vin,simOptions);
                case 'nmos'
                    c = this.cNMOS(deviceParams,Vin,simOptions);
                case 'vsrc'
                    c = zeros(size(Vin)); %voltage sources for now add no capacitance
                case 'inductor'
                    c = this.cL(deviceParams,Vin,simOptions);
                case 'capacitor'
                    c = this.cC(deviceParams,Vin,simOptions);
                case 'resistor'
                    c = this.cR(deviceParams,Vin,simOptions);
                otherwise
                    error(strcat(deviceType,' is not a supported device type'))
            end
        end
        
       
        %interface to evaluate PMOS device, must be overriden by specific
        %model class.
        function [i,cParam] = iPMOS(this,deviceParams,Vin,simOptions)
            i = [];
            cParam = [];
        end
            
        %interface to evaluate NMOS device, must be overriden by specific
        %model class.    
        function [i,cParam] = iNMOS(this,deviceParams,Vin,simOptions)
            i =[];
            cParam = [];
        end
        %interface to evaluate inductor device, must be overriden by specific
        %model class.        
        function [i,cParam] = iL(this,deviceParams,Vin,simOptions)
            i = [];
            cParam = [];
        end
        %interface to evaluate capacitor device, must be overriden by specific
        %model class.        
        function [i,cParam] = iC(this,deviceParams,Vin,simOptions)
            i = [];
            cParam = [];
        end
        %the default bahaviour is:         
        function [i,cParam] = iR(this,deviceParams,Vin,simOptions)
            i = Vin./deviceParams.R;
            cParam = 0;
        end
        
        %interface to evaluate PMOS device, must be overriden by specific
        %model class.
        function c = cPMOS(this,deviceParams,Vin,simOptions)
            c = [];
        end
            
        %interface to evaluate NMOS device, must be overriden by specific
        %model class.    
        function c = cNMOS(this,deviceParams,Vin,simOptions)
            c =[];
        end
        %interface to evaluate inductor device, must be overriden by specific
        %model class.        
        function c = cL(this,deviceParams,Vin,simOptions)
            c = [];
        end
        %interface to evaluate capacitor device, must be overriden by specific
        %model class.        
        function c = cC(this,deviceParams,Vin,simOptions)
            c = [];
        end
        %interface to evaluate resistor device, must be overriden by specific
        %model class.        
        function c = cR(this,deviceParams,Vin,simOptions)
            c = [];
        end
        
        %interface to evaluate residual current in NMOS device, must be
        %overriden by specific model class.
        function [ires,J] = iresNMOS(this,deviceParams,V,I,numCircuitNodes,simOptions)
            ires = [];
            J = [];
        end
        
        %interface to evaluate residual current in PMOS device, must be
        %overriden by specific model class.
        function [ires,J] = iresPMOS(this,deviceParams,V,I,numCircuitNodes,simOptions)
            ires = [];
            J = [];
        end
        
                %interface to evaluate inductor device, must be overriden by specific
        %model class.        
        function [ires,J] = iresL(this,deviceParams,V,I,numCircuitNodes,simOptions)
            ires = [];
            J = [];
        end
        %interface to evaluate capacitor device, must be overriden by specific
        %model class.        
        function [ires,J] = iresC(this,deviceParams,V,I,numCircuitNodes,simOptions)
            ires = [];
            J = [];
        end
        %the default bahaviour is:         
        function [ires,J] = iresR(this,deviceParams,V,I,numCircuitNodes,simOptions)
            ires = V./deviceParams.R - I;
            J = ires.dx;
            
        end
        
         
    end
end

