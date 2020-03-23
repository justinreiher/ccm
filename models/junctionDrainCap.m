function [ Cj] = junctionDrainCap( Vdb,params )
%junctionCap computes the junction capacitance between drain and body (substrate) 
%   INPUTS:
%   Vdb - is the Vd - Vb for the device
%
%   OUTPUTS:
%   Cj - is the total computed junction capacitance of the device


% Justin Reiher 2017/07/30


%from 45nm model card, celcius to kelvin
TNOM = params.TNOM; %degrees Kelvin

T = params.Tjun;
%bottom capacitance 
CJD = params.CjdTnom*(1+params.Tcj*(T-TNOM));
PBD = params.PbdTnom - params.Tpb*(T-TNOM);

%sidewall capacitance
CJSWD = params.CjswdTnom + params.Tcjsw*(T-TNOM);
PBSWD = params.PbswdTnom - params.Tpbsw*(T-TNOM);

Xj = params.Xj; %[cm] from model card

%sidewall to gate capacitance
CJSWGD = params.CjswgdTnom*(1+params.Tcjswg*(T-TNOM));
PBSWGD = params.PbswgdTnom - params.Tpbswg*(T-TNOM);

% Lactive = Lg + XL - 2*dL
Lactive = params.Lgdr +params.Xl -2*params.dLg;

% Wactive = W + XW - 2*dW
Wactive = params.W + params.Xw - 2*params.Wint;
Ad = Lactive.*Wactive;
Pd = Xj*(Wactive + 2*Lactive);

MJD = 0.5;
MJSWD = 0.33;
MJSWGD = 0.33;

indNeg = Vdb < 0;

%Cjb = Vdb*0;
%Cjswg = Vdb*0;
%Cjsw  = Vdb*0;

%if all of Vdb < 0 then do the following

if(all(indNeg == 1))
    %Computes the parts where Vdb < 0
    CjbN = CJD*(1-Vdb(indNeg)/PBD).^(-MJD).*Ad(indNeg)/1e4;
    CjswgN = CJSWGD*(1-Vdb(indNeg)/PBSWGD).^(-MJSWGD).*Wactive(indNeg)/1e2;
    CjswN =  CJSWD*(1-Vdb(indNeg)/PBSWD).^(-MJSWD).*Pd(indNeg)/1e4 ;
    Cj = CjbN + CjswgN + CjswN;
elseif(all(indNeg == 0))
    %Computes the parts where Vdb > 0
    CjbP = CJD*(1+MJD*Vdb(~indNeg)/PBD).*Ad(~indNeg)/1e4;
    CjswgP = CJSWGD*(1+MJSWGD*Vdb(~indNeg)/PBSWGD).*Wactive(~indNeg)/1e2;
    CjswP = CJSWD*(1+MJSWD*Vdb(~indNeg)/PBSWD).*Pd(~indNeg)/1e4;
    Cj =  CjbP + CjswgP + CjswP;
else
    %Computes the parts where Vdb < 0
    CjbN = CJD*(1-Vdb(indNeg)/PBD).^(-MJD).*Ad(indNeg)/1e4;
    CjswgN = CJSWGD*(1-Vdb(indNeg)/PBSWGD).^(-MJSWGD).*Wactive(indNeg)/1e2;
    CjswN =  CJSWD*(1-Vdb(indNeg)/PBSWD).^(-MJSWD).*Pd(indNeg)/1e4 ;
    CjN = CjbN + CjswgN + CjswN;
    %Computes the parts where Vdb > 0
    CjbP = CJD*(1+MJD*Vdb(~indNeg)/PBD).*Ad(~indNeg)/1e4;
    CjswgP = CJSWGD*(1+MJSWGD*Vdb(~indNeg)/PBSWGD).*Wactive(~indNeg)/1e2;
    CjswP = CJSWD*(1+MJSWD*Vdb(~indNeg)/PBSWD).*Pd(~indNeg)/1e4;
    CjP =  CjbP + CjswgP + CjswP;
    
    Cj = Vdb*0;
    
    
    Cj(indNeg) = CjN;
    Cj(~indNeg) = CjP;
    
end

end

