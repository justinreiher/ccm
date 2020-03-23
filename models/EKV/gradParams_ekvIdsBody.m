%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a simplified EKV model which computes the current flowing through
% a transistor device model with the following fit parameters:
%   - alpha: Parameter to model Drain Induced Barrier Lowering
%   - beta:  Paramater to model Drain Induced Barrier Lowering
%   - Vth:   Parameter that models the threshold voltage for the device
%   - W:     Parameter that models the width of the device
%   - Lg:    Parameter that models the length of the device
%   - Io:    The benchmark current
%
% Returns gradParams_Ids = [partial Ids/partial Io ; partial Ids / partial alpha ; partial
% Ids / partial beta ; partial Ids / partial Vth]
%
% Derived by Mark R. Greenstreet and Justin Reiher at University of British
% Columbia
%
% Version 1.0
% March 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gradParams_Ids = gradParams_ekvIds(params,biasVoltage)
    alpha = params.alpha;
    beta  = params.beta;
    Vth   = params.Vth;
    W     = params.W;
    Lg    = params.Lg;
    Io    = params.Io;
    gamma = params.gamma;
    phi   = params.phi;
    type  = params.type;
    
    Vd = biasVoltage(1,:);
    Vg = biasVoltage(2,:);
    Vs = biasVoltage(3,:);
    Vb = biasVoltage(4,:);
    
    phiM = sqrt(phi + type*(Vs-Vb)) - sqrt(phi);
    Vmod = gamma*(phiM);
    
    u = alpha*(Vg+beta*Vd-Vs-(Vth+Vmod));
    v = alpha*(Vg+beta*Vs-Vd-(Vth+Vmod));
    
    partialIo = W*(log(1+exp(u))-log(1+exp(v)));
    partialAlpha = W*Io*(exp(u).*(Vg+beta*Vd-Vs-Vth)./(1+exp(u)) - exp(v).*(Vg+beta*Vs-Vd-Vth)./(1+exp(v)));
    partialBeta = W*Io*(exp(u)*alpha.*Vd./(1+exp(u)) - exp(v)*alpha.*Vs./(1+exp(v)));
    partialVth = W*Io*(-alpha*exp(u)./(1+exp(u)) + alpha*exp(v)./(1+exp(v)));
    

    partialGamma = W*Io*((-alpha*exp(u).*phiM)./(1+exp(u)) + ...
                           (alpha*exp(v).*phiM)./(1+exp(v)));
    partialPhi   = W*Io*((-alpha*exp(u).*gamma./(1+exp(u))).*(1./(2*sqrt(phi+type*(Vs-Vb)))-1/(2*sqrt(phi))) + ...
                          (alpha*exp(v).*gamma./(1+exp(v))).*(1./(2*sqrt(phi+type*(Vs-Vb)))-1/(2*sqrt(phi))));
    
    gradParams_Ids = [partialIo;
                    partialAlpha;
                    partialBeta;
                    partialVth;
                    partialGamma;
                    partialPhi];
end
