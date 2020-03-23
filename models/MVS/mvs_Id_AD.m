function [Id_ret, intermediate_value,J] = mvs_Id_AD ( params, biasVoltages,IdspreComputed) %Vd,Vg,Vs,Vb )
% Symmetrical Short-channel MOSFET Compact Model (Revised)
%   for Automatic Differentiation

% This is a Revision of the Symmetrical Short-channel MOSFET
% Model originally written by Dimitri Antoniadis, MIT:
% https://nanohub.org/publications/15/4

% This revision is developed to be able to run automatic
% differentiation over the model. 
% 
% Currently there are two suspicious functions that are not 
% differentiable: if-then-else and abs.

% Yan Peng, UBC, 2016/09/28
% Justin Reiher, UBC, 2017/06/08 *modified to use just vectors, no AD
% Justin Reiher, UBC 2017/06/20 re-introduced AD, but still just returning
%                Id as a value

% Parameters
    version = params.version;
    type = params.type;
    W = params.W;
    Lgdr = params.Lgdr;
    dLg = params.dLg;
    Cg = params.Cg;
    etov = params.etov;
    delta = params.delta;
    n0 = params.n0;
    Rs0 = params.Rs0;
    Rd0 = params.Rd0;
    Cif = params.Cif;
    Cof = params.Cof;
    vxo = params.vxo*1e7;
    mu = params.mu;
    beta = params.beta;
    Tjun = params.Tjun;
    phib = params.phib;
    gamma = params.gamma;
    Vt0 = params.Vt0;
    alpha = params.alpha;
    mc = params.mc;
    CTM_select = params.CTM_select;
    CC = params.CC;
    nd = params.nd;
    zeta = params.zeta;

    % SMALL_VALUE
    SMALL_VALUE = 1e-10;
    % LARGE_VALUE
    LARGE_VALUE = 40;
    
    Rs = 1e-4./ W * Rs0;

    %bias voltages
    Vd = biasVoltages(1,:);
    Vg = biasVoltages(2,:);
    Vs = biasVoltages(3,:);
    Vb = biasVoltages(4,:);

    % Initial values
    %Vdi = Vd;
    %Vsi = Vs;
    dir = zeros(size(Vs));
    
    Vdsraw = type*(biasVoltages(1,:) - biasVoltages(3,:));
    
    s_bigger_than_d = find(Vdsraw(:) < 0);
    s_smaller_than_d = find(Vdsraw(:) >= 0);
    
    VsTemp = Vs;
    VdTemp = Vd;
    
    VsTemp(s_bigger_than_d) = Vd(s_bigger_than_d);
    VdTemp(s_bigger_than_d) = Vs(s_bigger_than_d);
    
    dir(s_smaller_than_d) = 1;
    dir(s_bigger_than_d) = -1;

    % Voltage definitions
    Vgsraw = type*( biasVoltages(2,:) - biasVoltages(3,:));
    Vgdraw = type*( biasVoltages(2,:) - biasVoltages(1,:));

    %Initialize the variables
%     Vds = Vd - Vs;
%     Vbs = Vb - Vs;
%     Vgs = Vg - Vs;
%     Vdsi = Vdi-Vsi;
%     Vbsi = Vb - Vsi;
%     Vgsi = Vgsraw;

%     dir = zeros(size(Vs.x));
% 
%     VsTemp = Vs;
%     VdTemp = Vd;
% 
%     s_bigger_than_d = find(Vgsraw(:) >= Vgdraw(:));
%     s_smaller_than_d = find(Vgsraw(:) < Vgdraw(:));
%     
%     VsTemp(s_smaller_than_d) = Vd(s_smaller_than_d);
%     VdTemp(s_smaller_than_d) = Vs(s_smaller_than_d);
%     
%     dir(s_bigger_than_d) = 1;
%     dir(s_smaller_than_d) = -1;
    
    Vds = type * (VdTemp - VsTemp);
    Vgs = type * (Vg - VsTemp);
    Vbs = type * (Vb - VsTemp);
    Vdsi = type * (Vd - VsTemp);
    Vgsi = Vg - VsTemp;
    Vbsi = type * (Vb - VsTemp);
    
    if(any(Vds == 0))
        Vds = Vds + 1e-100;
    end
    
%     the above replaces what below was supposed to do
% what happens in the code below is that if both s_bigger_than_d and
% s_smaller_than_d are non-empty, the first set of Vds gets assigned and
% then when it tries to find the indexes of s_smaller_than_d it is an out
% of bounds index error. But I can't initialize an empty variable of the
% right size because the way AD creates variables doesn't allow it to place
% the values as below
    
%     if ~isempty(s_bigger_than_d)
%         Vds(s_bigger_than_d)   = type * ( biasVoltages(1,(s_bigger_than_d)) - biasVoltages(3,(s_bigger_than_d)) );
%         Vgs(s_bigger_than_d)   = type * ( biasVoltages(2,(s_bigger_than_d)) - biasVoltages(3,(s_bigger_than_d)) );
%         Vbs(s_bigger_than_d)   = type * ( biasVoltages(4,(s_bigger_than_d)) - biasVoltages(3,(s_bigger_than_d)) );
%         Vdsi(s_bigger_than_d)    = type * ( biasVoltages(1,(s_bigger_than_d)) - biasVoltages(3,(s_bigger_than_d)) );
%         Vgsi(s_bigger_than_d)    = Vgsraw(s_bigger_than_d);
%         Vbsi(s_bigger_than_d)    = type * ( biasVoltages(4,(s_bigger_than_d)) - biasVoltages(3,(s_bigger_than_d)) );
%         dir(s_bigger_than_d)   = 1;
%     elseif ~isempty(s_smaller_than_d)
%         Vds(s_smaller_than_d)    = type * ( biasVoltages(3,(s_smaller_than_d)) - biasVoltages(1,(s_smaller_than_d)) );
%         Vgs(s_smaller_than_d)    = type * ( biasVoltages(2,(s_smaller_than_d)) - biasVoltages(1,(s_smaller_than_d)) );
%         Vbs(s_smaller_than_d)    = type * ( biasVoltages(4,(s_smaller_than_d)) - biasVoltages(1,(s_smaller_than_d)) );
%         Vdsi(s_smaller_than_d)   = type * ( biasVoltages(3,(s_smaller_than_d)) - biasVoltages(1,(s_smaller_than_d)) );
%         Vgsi(s_smaller_than_d)   = Vgdraw(s_smaller_than_d);
%         Vbsi(s_smaller_than_d)   = type * ( biasVoltages(4,(s_smaller_than_d)) - biasVoltages(1,(s_smaller_than_d)) );
%         dir(s_smaller_than_d)    = -1;
%     else 
%         biasVoltages
%         error('something went wrong...');
%     end


     %dvds =0;
     %dvg = 0;
%         fprintf('Iterative ')
  %  tic()
%  dvds = 0; dvg =0;
  % [dvds,dvg, Ids] = virtualSource(params, [Vds.x;Vgs.x;Vbs.x;Vdsi.x;Vs.x;Vd.x;Vg.x],type);
dvg = IdspreComputed.*Rs;
dvds = 2*dvg;
   IdsAD = IdspreComputed .*dir.*type;
  dvdsAD = dvds;
   dvgAD = dvg;
 % [dvdsAD,dvgAD] = virtualSourceAD(params,[Vds;Vgs;Vbs;Vdsi;Vs;Vd;Vg],dvds,dvg,type);
  % toc()
  % [dV,Id_ret] = wrapper(Vds,Vgs,Vbs,dvds,dvg,Vt0,gamma,phib,delta,aphit,phit,Cg,nphit,vxo,Leff,mu,beta,type,dir);
    Vdsi=Vds-dvdsAD;
    Vgsi=Vgs-dvgAD;
    Vbsi=Vbs-dvgAD;
  %  toc()
        

    % Parasitic element definition
%    Rs = 1e-4./ W * Rs0;                                      % s-terminal resistance [ohms]
%    Rd = Rs;                                                 % d-terminal resistance [ohms] For symmetric source and drain Rd = Rs.

    % d-terminal resistance [ohms] {Uncomment for asymmetric source and drain resistance.}
    % Rd      = 1e-4/ W * Rd0;

    Cofs = 0*( 0.345e-12/ etov ) * dLg/ 2.0 + Cof;  % s-terminal outer fringing cap [F/cm]
    Cofd = 0*( 0.345e-12/ etov ) * dLg/ 2.0 + Cof;  % d-terminal outer fringing cap [F/cm]
    Leff = Lgdr - dLg;                            % Effective channel length [cm]. After subtracting overlap lengths on s and d side

    kB = 8.617e-5;                                % Boltzmann constant [eV/K]
    phit = kB*Tjun;                               % Thermal voltage, kT/q [V]
    me = (9.1e-31) * mc;                          % Carrier mass [Kg]
    n = n0 + nd * Vdsi;                           % Total subthreshold swing factor taking punchthrough into account [unit-less]
    nphit = n * phit;                             % Product of n and phit [used as one variable]
    aphit = alpha * phit;                         % Product of alpha and phit [used as one variable]

    % Correct Vgsi and Vbsi
    % Vcorr is computed using external Vbs and Vgs but internal Vdsi, Qinv and Qinv_corr are computed with uncorrected Vgs, Vbs and corrected Vgs, Vbs respectively.
    Vtpcorr = Vt0 + gamma * (sqrt(abs(phib - Vbs))- sqrt(phib))- Vdsi * delta; % Calculated from extrinsic Vbs
    eVgpre = exp(( Vgs - Vtpcorr )/ ( aphit * 1.5 ));                          % Calculated from extrinsic Vgs
    FFpre = 1.0./ ( 1.0 + eVgpre );                                            % Only used to compute the correction factor
    ab = 2 * ( 1 - 0.99 * FFpre ) * phit;
    Vcorr = ( 1.0 + 2.0 * delta ) * ( ab./ 2.0 ) .* ( exp( -Vdsi./ ab ));         % Correction to intrinsic Vgs
    Vgscorr = Vgsi + Vcorr;                                             % Intrinsic Vgs corrected (to be used for charge and current computation)
    Vbscorr = Vbsi + Vcorr;                                             % Intrinsic Vgs corrected (to be used for charge and current computation)
    Vt0bs = Vt0 + gamma .* (sqrt( abs( phib - Vbscorr)) - sqrt( phib )); % Computed from corrected intrinsic Vbs
    Vt0bs0 = Vt0 + gamma .* (sqrt( abs( phib - Vbsi)) - sqrt( phib ));   % Computed from uncorrected intrinsic Vbs
    Vtp = Vt0bs - Vdsi .* delta - 0.5 * aphit;                           % Computed from corrected intrinsic Vbs and intrinsic Vds
    Vtp0 = Vt0bs0 - Vdsi .* delta - 0.5 * aphit;                         % Computed from uncorrected intrinsic Vbs and intrinsic Vds
    eVg = exp(( Vgscorr - Vtp )/ ( aphit ));                            % Compute eVg factor from corrected intrinsic Vgs
    FF = 1.0./ ( 1.0 + eVg );
    eVg0 = exp(( Vgsi - Vtp0 )/ ( aphit ));                             % Compute eVg factor from uncorrected intrinsic Vgs
    FF0 = 1.0./ ( 1.0 + eVg0 );


    Qref = Cg * nphit;                                                  % F/cm^2*V=Coulombs/cm^2,   F = Coulombs/V
    eta = ( Vgscorr - ( Vt0bs - Vdsi .* delta - FF * aphit ))./ ( nphit ); % Compute eta factor from corrected intrinsic Vgs and intrinsic Vds
    eta0 = ( Vgsi - ( Vt0bs0 - Vdsi .* delta - FF * aphit ))./ ( nphit );  % Compute eta0 factor from uncorrected intrinsic Vgs and internal Vds.
                                                                         % Using FF instead of FF0 in eta0 gives smoother capacitances.

    % %Charge at VS in saturation (Qinv)
    Qinv_corr = Qref.*log_one_plus_exp(eta, 0);
    Qinv = Qref.*log_one_plus_exp(eta0, 0);
    
    
    %Transport equations
    vx0 = vxo;
    Vdsats = vx0 .* Leff./ mu;
    Vdsat = Vdsats .* ( 1.0 - FF ) + phit .* FF;   % Saturation drain voltage for current    
    Vdratio = abs( Vdsi)./ Vdsat;
    
    Vdbeta = power( Vdratio, beta);
    Vdbetabeta = power( 1.0 + Vdbeta, 1.0/ beta);
    Fsat = Vdratio ./ Vdbetabeta;                 % Transition function from linear to saturation.
                                                 % Fsat = 1 when Vds>>Vdsat; Fsat= Vds when Vds<<Vdsat
    %Total drain current
    Id = Qinv_corr .* vx0 .* Fsat .* W.*type.*dir;  % Coulombs/s = A
    
    
    Id_ret = Id;
   % dI = Id.dx;
    
%     deltaVs = abs(Id.x*Rs) - abs(Vs.x);
%     deltaVd = abs(Id.x*Rd) - abs(Vd.x);
%     
   

    

    

    %%% Total capacitance, on average, each transistor should provide 2.
    % First try:
    % C = W*Lgdr*Cg*4;

    % Second try:
    % Cd = (type == 1)*0.96e-16 + (type == -1)*1.34e-16;
    % Cs = (type == 1)*0.96e-16 + (type == -1)*1.34e-16;
    % Cg = (type == 1)*0.96e-16 + (type == -1)*1.34e-16;
    % C = [Cd; Cs; Cg];

    % want to solve:
    % [VsInt] =  [1-Rs*d/dVs(Id)   -Rs*d/dVd(Id) ]  [dVs/dI]  
    % [VdInt] =  [Rd*d/dVs(Id)     1+Rd*d/dVd(Id)]  [dVd/dI]
    
    % dV is not used at this moment.
dV=Id;        
       
     Vsint=Vs+Id.*(Rs0*1e-4./W).*dir;
     Vdint=Vd-Id.*(Rd0*1e-4./W).*dir;
     Vgsraw=type*(Vg-Vsint);
     Vgdraw=type*(Vg-Vdint);
    
    

    % Third try
    intermediate_value = struct('Vb', biasVoltages(4), 'LARGE_VALUE', LARGE_VALUE, 'SMALL_VALUE', SMALL_VALUE, 'Qinv', Qinv, 'eta0', eta0, 'phit', phit, 'Vsi', Vsint, 'Vdi', Vdint, 'Vbsi', Vbsi, 'Vdsi', Vdsi, 'Vgsi', Vgsi, 'Vgsraw', Vgsraw, 'Vgdraw', Vgdraw, 'FF0', FF0, 'aphit', aphit, 'nphit', nphit, 'me', me, 'vx0', vx0, 'Leff', Leff, 'Vt0bs0', Vt0bs0, 'FF', FF, 'Qref', Qref, 'dir', dir, 'Cofs', Cofs, 'Cofd', Cofd, 'Fsat', Fsat);
end

