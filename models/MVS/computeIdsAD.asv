function IdsRes = computeIdsAD( params,mvsAD)
%virtualSource Computes via an iterative solver technique the voltages for
%dvds and dvg used to accurately compute Ids of the device. This code is
%vectorized so that it can compute these voltages in parallel for n devices
%of the same type (either nmos or pmos).

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
    
    kB = 8.617e-5;                                % Boltzmann constant [eV/K]
    phit = kB*Tjun;                               % Thermal voltage, kT/q [V]
    me = (9.1e-31) * mc;                          % Carrier mass [Kg]
    
    
    Cofs = 0*( 0.345e-12/ etov ) * dLg/ 2.0 + Cof;  % s-terminal outer fringing cap [F/cm]
    Cofd = 0*( 0.345e-12/ etov ) * dLg/ 2.0 + Cof;  % d-terminal outer fringing cap [F/cm]
    Leff = Lgdr - dLg;  % Effective channel length [cm]. After subtracting overlap lengths on s and d side
    
    Vd = mvsAD(1,:);
    Vg = mvsAD(2,:);
    Vs = mvsAD(3,:);
    Vb = mvsAD(4,:);
    IdsAD = mvsAD(5,:);
       
    Vdsraw = type*(mvsAD(1,:) - mvsAD(3,:));
    
    s_bigger_than_d = find(Vdsraw(:) < 0);
    s_smaller_than_d = find(Vdsraw(:) >= 0);
    
    VsTemp = Vs;
    VdTemp = Vd;
    
    VsTemp(s_bigger_than_d) = Vd(s_bigger_than_d);
    VdTemp(s_bigger_than_d) = Vs(s_bigger_than_d);
    
    dir(s_smaller_than_d) = 1;
    dir(s_bigger_than_d) = -1;

    % Voltage definitions
    Vgsraw = type*( mvsAD(2,:) - mvsAD(3,:));
    Vgdraw = type*( mvsAD(2,:) - mvsAD(1,:));

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
    
    if(any(Vds == 0))
        Vds = Vds + 1e-100;
    end
    
    Rs = 1e-4./ W * Rs0;
    Rd = Rs;

    dvg=IdsAD.*Rs;
    dvd=IdsAD.*Rd;
    dvds = dvg + dvd;
    
    Vdsi=Vds-dvds;
    Vgsi=Vgs-dvg;
    Vbsi=Vbs-dvg;
    
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
    IdsRes = Qinv_corr .* vx0 .* Fsat .* W.*type.*dir - IdsAD;  % Coulombs/s = A
    
    
    
    
  
    

end

