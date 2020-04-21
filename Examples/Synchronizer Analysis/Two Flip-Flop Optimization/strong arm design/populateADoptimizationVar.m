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

% This version of the function is to optimize a 2-flip flop strongARM latch
% synchronizer where the coupling inverters are free to be individually
% sized.

% buffer = 1,2 strongARM latch = 2 to 18
dev_indN = [1,2, 3,4, 5,6,7,8, 9,10,11,12];
dev_indP = [1,2, 3,4, 5,6, 7,8,9,10];

numVar = length(dev_indN) + length(dev_indP);

v = 1:length(voltageState);
nBuffer1 = 1+length(v):1+length(v);
nBuffer2 = 1+length(v)+length(nBuffer1):1+length(v)+length(nBuffer1);
nFF = 1+length(v)+length(nBuffer1)+length(nBuffer2):length(v)+length(dev_indN);
pBuffer1 = 1+length(v)+length(dev_indN):1+length(v)+length(dev_indN);
pBuffer2 = 1+length(v)+length(dev_indN)+length(pBuffer1):1+length(v)+length(dev_indN)+length(pBuffer1);
pFF = 1+length(v)+length(dev_indN)+length(pBuffer1)+length(pBuffer2):length(v)+length(dev_indN)+length(dev_indP);


AD_f = hessianinit([voltageState;nDevices(dev_indN)';pDevices(dev_indP)']);

adVar = [AD_f(v);
         AD_f(nBuffer1);
         AD_f(nBuffer2);
         AD_f(nFF);
         AD_f(nFF);
         AD_f(pBuffer1);
         AD_f(pBuffer2);
         AD_f(pFF);
         AD_f(pFF);];




end

