function stop = nestedBisectionStop(tbOptions,numStatesPerCCT,Vin,tCrit,t)
%nestedBisectionStop is a function which determines the stopping condition upon
%which to terminate the nested bisection algorithm. For example when the
%before last stage of the synchronizer has resolved metastability by
%tCrit. When this function gets called the simulation time has exceeded the
%user specified tCrit time and this function evaluates the selected
%criteria to return TRUE if the nested bisection alogirthm is to halt,
%otherwise returns FALSE.
%
%Inputs: numParallelCCTs - the number of parallel circuits being simulated.
%        numStatesPerCCT - the number of states per circuit
%        Vin             - the integrated voltage state vector
%        tCrit           - the critical time
%        t               - the sim time

numParallelCCTs = tbOptions.numParallelCCTs;
railLimit = tbOptions.digitalSimOptions.threshold;
vdd = tbOptions.vdd;
failHigh = (1-railLimit)*vdd;
failLow  = railLimit*vdd;

metaResolveNode = 14-4;
outputNode = 7-4;

stop = false;
hCond= false;
lCond = false;
ind = find(t>= tCrit);
ind = ind(1);

for i = 1:numParallelCCTs
    ViQmeta = Vin(ind,metaResolveNode + (i-1)*numStatesPerCCT);
    ViQEndmeta = Vin(end,metaResolveNode + (i-1)*numStatesPerCCT);
    ViQEnd     = Vin(end,outputNode + (i-1)*numStatesPerCCT);
    if( (ViQmeta < failHigh) && (ViQEndmeta > failHigh) && (ViQEnd > failHigh))
        hCond = true;
    end
    if( (ViQmeta > failLow) && (ViQEndmeta < failLow) && (ViQEnd < failLow))
        lCond = true;
    end
end

if( hCond && lCond)
    stop = true;
end
end

