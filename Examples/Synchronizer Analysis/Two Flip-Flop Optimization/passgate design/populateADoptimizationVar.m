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

% This version of the function is to optimize a 2-flip flop passgate
% synchronizer where the coupling inverters are free to be individually
% sized.

% buffer = 1, passgate latch = 2 to 6
dev_indN = [1,2,3,4,5,6];
dev_indP = [1,2,3,4,5,6];

numVar = length(dev_indN) + length(dev_indP);

vCol = 1:length(voltageState);
nColBuf = length(voltageState)+1:length(voltageState)+1;
nCol    = length(voltageState)+2:length(voltageState)+length(dev_indN);

pColBuf = length(voltageState)+length(dev_indN) + 1: length(voltageState)+length(dev_indN) + 1;
pCol    = length(voltageState)+length(dev_indN) + 2: length(voltageState)+length(dev_indN) + length(dev_indP);

offsetV = length(voltageState);
bufOffset = 1;
offsetN = length(dev_indN)-1;
offsetP = length(dev_indP)-1;

AD_f = hessianinit([voltageState;nDevices(dev_indN)';pDevices(dev_indP)']);

vState = AD_f(1:length(voltageState));
buffN   = AD_f(offsetV+1);
latchN  = AD_f(offsetV+2:offsetN);
buffP   = AD_f(offsetN + 1);
latchP  = AD_f(offsetN+2:end);

syncN = [buffN;latchN;latchN;latchN;latchN];
syncP = [buffP;latchP;latchP;latchP;latchP];

syncState = diag(ones(1,length(voltageState)));

bufferN = 1;
latchN = diag(ones(1,length(dev_indN)-1));

bufferP = 1;
latchP = diag(ones(1,length(dev_indP)-1));

map = zeros(length(voltageState)+length(nDevices)+length(pDevices),length(voltageState)+length(dev_indN) + length(dev_indP));

%voltage State
map(1:offsetV,vCol) = syncState;
%n-devices
%buffer
map(1+offsetV:offsetV+bufOffset,nColBuf) = bufferN;
%latches
map(1+offsetV+bufOffset:offsetV+bufOffset+offsetN,nCol) = latchN;
map(1+offsetV+bufOffset+offsetN:offsetV+bufOffset+2*offsetN,nCol) = latchN;
map(1+offsetV+bufOffset+2*offsetN:offsetV+bufOffset+3*offsetN,nCol) = latchN;
map(1+offsetV+bufOffset+3*offsetN:offsetV+bufOffset+4*offsetN,nCol) = latchN;
%p-devices
map(1+offsetV+bufOffset+4*offsetN:offsetV+bufOffset+4*offsetN+bufOffset,pColBuf) = bufferP;
map(1+offsetV+bufOffset+4*offsetN+bufOffset:offsetV+bufOffset+4*offsetN+bufOffset+offsetP,pCol) = latchP;
map(1+offsetV+bufOffset+4*offsetN+bufOffset+offsetP:offsetV+bufOffset+4*offsetN+bufOffset+2*offsetP,pCol) = latchP;
map(1+offsetV+bufOffset+4*offsetN+bufOffset+2*offsetP:offsetV+bufOffset+4*offsetN+bufOffset+3*offsetP,pCol) = latchP;
map(1+offsetV+bufOffset+4*offsetN+bufOffset+3*offsetP:end,pCol) = latchP;

adVar = map*AD_f;





end

