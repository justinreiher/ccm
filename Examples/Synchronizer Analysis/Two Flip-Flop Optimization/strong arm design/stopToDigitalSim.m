 function [value,isterminal,direction] = stopToDigitalSim(t,Vin,vDot,digitalSimOptions)
        
        outputs = [7,8,9,10,11,12,13,14] - 4;
 
        mask = zeros(1,20);
        mask(outputs) = 1;
        nodesPerCCT = length(mask);
        numParallelCCTs = length(Vin)/nodesPerCCT;
        
        mask = repmat(mask,1,numParallelCCTs);
        condMask = sum(mask);
        mask = logical(mask);
        
        vdd = digitalSimOptions.vdd;
        minSimTime = digitalSimOptions.minSimTime;
        threshold = digitalSimOptions.threshold;
        stopGain = digitalSimOptions.stopGain;
        
        stopCondition = sum(tanh(stopGain*(vdd*(1-threshold)*ones(1,condMask)*1.05-vdd/2).^2 ...
          -(threshold - vdd/2)^2));
 
        if(t == minSimTime)
            t = minSimTime + 1e-6;
        end
        
        a = stopGain*sign(t-minSimTime);
        isterminal = 1;
        direction = 1;
        value = sum(tanh(a*(Vin(mask)-vdd/2).^2-...
            (threshold-vdd/2)^2))-stopCondition;
        if value > 0 
            value = value + 0.5 - norm(vDot);
        end
        
        
    end