%%%
% Function for returning node charges.
% Yan Peng, UBC, 2016/10/24
%
% Justin REiher, UBC, 2020/03/24
% Modified function to return Capacitance - Requires AD

function [C] = mvs_c_AD(params, init, iv)

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

% Intermediate values
LARGE_VALUE = iv.LARGE_VALUE;
SMALL_VALUE = iv.SMALL_VALUE;
Qinv = iv.Qinv;
eta0 = iv.eta0;
phit = iv.phit;
Vsi = iv.Vsi;
Vdi = iv.Vdi;
Vbsi = iv.Vbsi;
Vdsi = iv.Vdsi;
Vgsi = iv.Vgsi;
Vgsraw = iv.Vgsraw;
Vgdraw = iv.Vgdraw;
FF0 = iv.FF0;
aphit = iv.aphit;
nphit = iv.nphit;
me = iv.me;
vx0 = iv.vx0;
Leff = iv.Leff;
Vt0bs0 = iv.Vt0bs0;
FF = iv.FF;
Qref = iv.Qref;
dir = iv.dir;
Cofs = iv.Cofs;
Cofd = iv.Cofd;
Fsat = iv.Fsat;

Vd = init(1,:);
Vg = init(2,:);
Vs = init(3,:);
Vb = init(4,:);



% Calculation of intrinsic charge partitioning factors (qs and qd)
Vgt	=	Qinv/ Cg;               % Use charge computed from uncorrected intrinsic Vgs

% Approximate solution for psis is weak inversion
if (gamma == 0)
    a = 1.0;
    % if (eta0 <= LARGE_VALUE)
    %     psis = phib + phit * ( 1.0 + log( log( 1.0 + SMALL_VALUE + exp( eta0 ))));
    % else
    %     psis = phib + phit * ( 1.0 + log( eta0 ));
    % end
    psis = phib + phit * ( 1.0 + log( log_one_plus_exp(eta0, SMALL_VALUE)));
else
    % if (eta0 <= LARGE_VALUE)
    %     psis = phib + ( 1.0 - gamma )/ ( 1.0 + gamma ) * phit * ( 1.0 + log( log( 1.0 + SMALL_VALUE + exp( eta0 ))));
    % else
    %     psis = phib + ( 1.0 - gamma )/ ( 1.0 + gamma ) * phit * ( 1.0 + log( eta0 ));
    % end
    psis = phib + ( 1.0 - gamma )/ ( 1.0 + gamma ) * phit * ( 1.0 + log( log_one_plus_exp(eta0, SMALL_VALUE)));
    a = 1.0 + gamma./ ( 2.0 * sqrt( abs( psis - ( Vbsi ))));
end
Vgta = Vgt./ a;               %   Vdsat in strong inversion
Vdsatq = sqrt( FF0 .* aphit .* aphit + Vgta .* Vgta);     %   Vdsat approx. to extend to weak inversion;
                                                       %   The multiplier of phit has strong effect on Cgd discontinuity at Vd=0.

% Modified Fsat for calculation of charge partitioning
% DD-NVSAT charge
Fsatq = abs( Vdsi./ Vdsatq )./ ( power( 1.0 + power( abs( Vdsi./ Vdsatq ), beta ), 1.0/ beta ));
x = 1.0 - Fsatq;


% MVS 1.0.1 release uses Tsividis model for implementing DD-NVSAT charges.
%  In MVS 1.1.0, we use the formulation from L. Wei's 2012 TED paper referenced in the documentation.
% Both models are identical.

% *** From Tsividis (DD-NVSAT)
% den     = 15 * ( 1 + x ) * ( 1 + x );
% qsc       = Qinv *(6 + 12 * x + 8 * x * x + 4 * x * x * x)/ den;
% qdc     = Qinv *(4 + 8 * x + 12 * x * x + 6 * x * x * x)/ den;
% qi      = qsc + qdc;              %   Charge in the channel

% *** From L. Wei (DD-NVSAT)
A = power((1-x),2)./(12.0*(1-(1-x)/2.0));
B = (5.0-2.0*(1-x))./(10.0-(1-x)/2.0);
qsc = Qinv.*(0.5-(1-x)/6.0+A.*(1-B));
qdc = Qinv.*(0.5-(1-x)/3.0+A.*B);
qi = qsc+qdc;

% QB charge
% kq       =   0.0;
tol = ( SMALL_VALUE * vxo/ 100.0 ) * ( SMALL_VALUE * vxo/ 100.0 ) * me/ ( 2 * constants('P_Q') );
if (Vdsi <= tol)
    kq2 = ( 2.0 * constants('P_Q')/ me * Vdsi * zeta)/ ( vx0 * vx0 ) * 10000.0;
    kq4 = kq2 .* kq2;

    % Uncomment following lines for linear potential profile
    qsb = Qinv .* (0.5 - kq2/12.0 +kq4/32.0);
    qdb = Qinv .* (0.5 - kq2/6.0 +3.0*kq4/32.0);

    % Uncomment following lines for parabolic potential profile
    %qsb    = Qinv * ( 0.5 - kq2/ 24.0 + kq4/ 80.0 );
    %qdb    = Qinv * ( 0.5 - 0.125 * kq2 + kq4/ 16.0 );
else
    kq = sqrt( 2.0 * constants('P_Q')/ me * Vdsi * zeta )/ vx0 * 100.0;
    kq2 = kq .* kq;

    % Uncomment following lines for linear potential profile
    qsb = Qinv .* (4.0*(kq2+1.0).*sqrt(kq2+1.0)-(6.0*kq2+4.0))./(3.0*kq2.*kq2);
    qdb = Qinv .* (2.0*(kq2-2).*sqrt(kq2+1.0)+4.0)./(3.0*kq2.*kq2);

    % Uncomment following lines for parabolic potential profile
     %qsb   = Qinv * ( asinh( kq )/ kq - ( sqrt( kq2 + 1.0 ) - 1.0 )/ kq2);
     %qdb   = Qinv * (( sqrt( kq2 + 1.0 )- 1.0 )/ kq2);
end


% Flag for classic or ballistic charge partitioning:
%   Ballistic blended with classic DD-NVSAT
%   Calculation of "ballistic" channel charge partitioning factors, qsb and qdb.
%   Here it is assumed that the potential increases parabolically from the
%   virtual source point, where Qinv_corr is known to Vds-dvd at the drain.
%   Hence carrier velocity increases linearly by kq (below) depending on the
%   efecive ballistic mass of the carriers.
if (CTM_select == 1)
    qs = qsc;
    qd = qdc;
else
    qs = qsc .* ( 1 - Fsatq .* Fsatq ) + qsb .* Fsatq .* Fsatq;
    qd = qdc .* ( 1 - Fsatq .* Fsatq ) + qdb .* Fsatq .* Fsatq;
end


% Body charge based on approximate surface potential (psis) calculation with delta=0 using psis=phib in Qb gives continuous Cgs, Cgd, Cdd in SI, while Cdd is smooth anyway.
Qb = -type * W .* Leff * ( Cg * gamma * sqrt( abs( psis - Vbsi )) + ( a - 1.0 )./ ( 1.0 * a ) .* Qinv .* ( 1.0 - qi ));

% DIBL effect on drain charge calculation.
% Calculate dQinv at virtual source due to DIBL only. Then:Correct the qd factor to reflect this channel charge change due to Vd
% Vt0bs0 and FF=FF0 causes least discontinuity in Cgs and Cgd but produces a spike in Cdd at Vds=0 (in weak inversion.  But bad in strong inversion)
etai = ( Vgsi - ( Vt0bs0 - FF * aphit ))./ ( nphit );
% if (etai <= LARGE_VALUE)
%     Qinvi = Qref * log( 1.0 + exp( etai ));
% else
%     Qinvi = Qref * etai;
% end
Qinvi = Qref.*log_one_plus_exp(etai, 0);
dQinv = Qinv - Qinvi;
dibl_corr = ( 1.0 - FF0 ) .* ( 1.0 - Fsatq ) .* qi .* dQinv./Qinv;
qd = qd - dibl_corr;


% Inversion charge partitioning to terminals s and d
Qinvs = type * Leff * (( 1 + dir ) .* qs + ( 1 - dir ) .* qd)/ 2.0;
Qinvd = type * Leff * (( 1 - dir ) .* qs + ( 1 + dir ) .* qd)/ 2.0;


% Outer fringing capacitance
Qsov = Cofs * ( init(2) - Vsi );
Qdov = Cofd * ( init(2) - Vdi );


% Inner fringing capacitance
Vt0x = Vt0 + gamma * ( sqrt( abs( phib - type * ( init(4) - Vsi ))) - sqrt(phib));
Vt0y = Vt0 + gamma * ( sqrt( abs( phib - type * ( init(4) - Vdi ))) - sqrt(phib));
Fs_arg = ( Vgsraw - ( Vt0x - Vdsi .* delta .* Fsat ) + aphit * 0.5 )./ ( 1.1 * nphit );
% if (Fs_arg <= LARGE_VALUE)
%     Fs = 1.0 + exp( Fs_arg );
%     FFx = Vgsraw - nphit * log( Fs );
% else
%     Fs = 0.0;                % Not used
%     FFx = Vgsraw - nphit * Fs_arg;
% end
FFx = Vgsraw - nphit .* log_one_plus_exp(Fs_arg, 0);
Fd_arg = ( Vgdraw - ( Vt0y - Vdsi .* delta .* Fsat ) + aphit * 0.5 )./ ( 1.1 * nphit );
% if (Fd_arg <= LARGE_VALUE)
%     Fd = 1.0 + exp( Fd_arg );
%     FFy = Vgdraw - nphit * log( Fd );
% else
%     Fd = 0.0;                % Not used
%     FFy = Vgdraw - nphit * Fd_arg;
% end
FFy = Vgdraw - nphit .* log_one_plus_exp(Fd_arg, 0);
Qsif = ( type * Cif + CC * Vgsraw ) .* FFx;
Qdif = ( type * Cif + CC * Vgdraw ) .* FFy;

% Partitioned charge
Qs = -W .* ( Qinvs + Qsov + Qsif );           %   s-terminal charge
Qd = -W .* ( Qinvd + Qdov + Qdif );           %   d-terminal charge
Qg = -( Qs + Qd + Qb );            %   g-terminal charge

Q = [Qd; Qg; Qs; Qb];




end