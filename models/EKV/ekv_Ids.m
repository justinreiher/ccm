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
% Derived by Mark R. Greenstreet and Justin Reiher at University of British
% Columbia
%
% Version 1.0
% March 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ids,i_f,i_r,n] = ekv_Ids(params,biasVoltage)
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
    
    VthMod = gamma*(sqrt(phi + type*(Vs-Vb)) - sqrt(phi));
    
    i_f = log_one_plus_exp(alpha*(Vg + beta*Vd-Vs-(Vth+VthMod)),0); %forward  normalized current
    i_r = log_one_plus_exp(alpha*(Vg + beta*Vs-Vd-(Vth+VthMod)),0); %backward normalized current
    %drain source current
    Ids = Io*W.*(i_f - i_r);
    %n is the slope factor, used for the capacitance
    n = 1.3;
    %n = (1 - gamma./(2*sqrt(Vg - Vth + (gamma/2 + sqrt(phi))^2))).^(-1);

end
