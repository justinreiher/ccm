function [adVar,numVar] = populateADoptimizationVar(voltageState,nDevices,pDevices)
%populateADoptimizationVar populates the AD variables for the computation
%of optimizing synchronizer circuits. This function is used to let the user
%decide what transistor widths are available to be decoupled for
%optimization purposes.
%   Inputs: voltageState  - the voltage state vector of the synchronizer
%           nDevices      - the vector of NMOS device widths in the design
%           pDevices      - the vector of PMOS device widths in the design
%   Outputs: adVAR is the AD variables which are used to compute the
%   gradient with respect to the desired transistor widths in the form of a
%   hessian object [voltageState;nDevices;pDevices] as a column vector


AD_f = hessianinit([voltageState;nDevices';pDevices']);
numVar = length(nDevices)+length(pDevices);

adVar = eye(length(AD_f))*AD_f;





end

