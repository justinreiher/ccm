% Template function for stop to digital simulation for the nested bisection algorithm.
% Implements the event function for MATLAB's ode solvers.
% Available paramters are user defined (or default) digitalSimOptions
% struct
% the time point t, the voltage state of the synchronizer at time t and the
% state derivative of the synchronizer at time t. It is upto the user to
% determine what information they need to determine a stopping condition
% and how they want to encorporate that into the ode solver.

function [value,isterminal,direction] = stopToDigitalSim(t,Vin,vDot,digitalSimOptions)

    value = -1;
    isterminal = 0;
    direction = 1;

end