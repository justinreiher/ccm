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
function stop = nestedBisectionStop(tbOptions,numStatesPerCCT,Vin,tCrit,t)

stop = true;

end

