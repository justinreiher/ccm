 function [value,isterminal,direction] = stopToDigitalSim(t,Vin,vDot)
        numOfNodes = length(Vin);
        vdd = 1.0;
        minSimTime = 6e-10;
        threshold = 0.1;
        stopGain = 30;
        
        stopCondition = sum(tanh(stopGain*(vdd*(1-threshold)*ones(1,numOfNodes)*1.05-vdd/2).^2 ...
          -(threshold - vdd/2)^2));
 
        if(t == minSimTime)
            t = minSimTime + 1e-6;
        end
        
        a = stopGain*sign(t-minSimTime);
        isterminal = 1;
        direction = 1;
        value = sum(tanh(a*(Vin-vdd/2).^2-...
            (threshold-vdd/2)^2))-stopCondition;
    end