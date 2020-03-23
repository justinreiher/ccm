function [Idlog,Id,Qs,Qd,Qg,Qb,Vdsi_out]=mvs_si_1_1_0(input_parms,coeff,bias_data)
% Symmetrical Short-Channel MOSFET model (VERSION=1.1.0)

% Returns the log of drain current, Id [A] and partitioned charges
% This model is only valid for Vg >~ Vg(psis=phif) where psis is the surface
% potential.  I.e range of validity is from onset of weak inversion through
% strong inversion.

% Original Dimitri Antoniadis, MIT, 09/17/10
% Modified, DAA 10/20/12
% Modified, DAA 07/01/13
% Modified SR 07/24/13 
% Modified SR 09/19/13
% Modified DAA 05/19/15
% Modified SR 07/18/15

version = 1.10; % version number

%fitted coefficients
Rs0=coeff(1);       % Access region resistance for s terminal [Ohms-micron]
Rd0=Rs0;            % Access region resistance for d terminal [Ohms-micron] {Generally Rs0=Rd0 for symmetric source and drain}
delta=coeff(3);     % Drain induced barrier lowering (DIBL) [V/V]
n0=coeff(4);        % Subthreshold swing factor [unit-less] {typically between 1.0 and 2.0}
nd = coeff(5);      % Punch-through factor [1/V]
vxo = coeff(6)*1e7; % Virtual-source injection velocity [cm/s]
mu = coeff(7);      % low field mobility [cm^2/Vs]
Vt0 = coeff(8);     % Threshold voltage [V]

%% input parameters known and not fitted.
type=input_parms(1);    % type of transistor. nFET type=1; pFET type=-1
W=input_parms(2);       % Transistor width [cm]
Lgdr=input_parms(3);    % Physical gate length [cm]. This is the designed gate length for litho printing.
dLg = input_parms(4);   % Overlap length including both source and drain sides [cm].
gamma=input_parms(5);   % Body-factor [sqrt(V)]
phib = input_parms(6);  % ~2*phif [V]
Cg = input_parms(7);    % Gate-to-channel areal capacitance at the virtual source [F/cm^2]
Cif = input_parms(8);   % Inner-fringing capacitance [F/cm]
Cof = input_parms(9);   % Outer-fringing capacitance [F/cm]
etov = input_parms(10); % Equivalent thickness of dielectric at S/D-G overlap [cm]
mc=input_parms(11);     % Effective mass of carriers relative to m0 [unitless]
Tjun=input_parms(12);   % Junction temperature [K].
beta=input_parms(13);   % Saturation factor. Typ. nFET=1.8, pFET=1.6
alpha=input_parms(14);  % Empirical parameter associated with threshold voltage shift between strong and weak inversion.
CTM_select=input_parms(15); % Parameter to select charge-transport model 
                            % if CTM_select = 1, then classic DD-NVSAT
                            % model is used; for CTM_select other than
                            % 1,blended DD-NVSAT and ballistic charge
                            % transport model is used.
                            
zeta = input_parms(16);     % Energy-transfer factor that lies between zero and unity. 
CC = 0*3e-13 ;              % Fitting parameter to adjust Vg-dependent inner fringe capacitances {not used in this version.}

me=9.1e-31*mc;          % Effective mass [Kg] invoked for ballistic charges
qe=1.602e-19;           % Elementary charge [Col.]
kB = 8.617e-5;          % Boltzmann constant [eV/K]
Cofs=0*(0.345e-12/etov)*dLg/2 + Cof;  % s-terminal outer fringing cap [F/cm]
Cofd=0*(0.345e-12/etov)*dLg/2 + Cof;  % d-terminal outer fringing cap [F/cm]
Leff = Lgdr-dLg;                    % Effective channel length [cm]

%%%%% ####### model file begins
% Direction of current flow:
% dir=+1 when "x" terminal is the source
% dir=-1 when "y" terminal is the source |

%% bias values
Vd_pre= bias_data(:,1);
Vg_pre= bias_data(:,2);
Vb_pre=bias_data(:,3);
Vs_pre=bias_data(:,4);
phit = kB*Tjun;

for len_bias=1:length(Vd_pre)
    Vd=Vd_pre(len_bias);
    Vg=Vg_pre(len_bias);
    Vb=Vb_pre(len_bias);
    Vs=Vs_pre(len_bias);
    dir=type*sign(Vd-Vs);
    
    Vds=abs(Vd-Vs);
    Vgs=max(type*(Vg-Vs),type*(Vg-Vd));
    Vbs=max(type*(Vb-Vs),type*(Vb-Vd));
    
    Vt0bs=Vt0+gamma*(sqrt(abs(phib-Vbs))-sqrt(phib));
    
    % Denormalize access resistances and allocate them the "source" and
    % drain according to current flow
    Rs=1e-4/W*(Rs0*(1+dir)+Rd0*(1-dir))/2;
    Rd=1e-4/W*(Rd0*(1+dir)+Rs0*(1-dir))/2;
    
    n=n0+nd*Vds;
    aphit = alpha*phit;
    nphit = n*phit;
    Qref=Cg*nphit;
    
    %%% Initial values for current calculation %%%%%%%%%%%%%%%%%%%%
    FF=1./(1+exp((Vgs-(Vt0bs-Vds.*delta-0.5*aphit))/(aphit)));
    Qinv_corr = Qref.*log(1+exp((Vgs-(Vt0bs-Vds.*delta-FF*aphit))./(nphit)));
    Qinv = Qref.*log(1+exp((Vgs-Vt0bs)./(nphit)));
    Rt=Rs + Rd + (Lgdr-dLg)./(W*Qinv*mu);
    vx0=vxo;
    Vdsats=W*Qinv.*vx0.*Rt;
    Vdsat=Vdsats.*(1-exp(-Qinv./Qref))+phit*exp(-Qinv./Qref);
    Fsat=(1-exp(-2*Vds./Vdsat))./(1+exp(-2*Vds./Vdsat));
    Idx = W*Fsat.*Qinv_corr.*vx0;
    Idxx=1e-15;
    dvg=Idx.*Rs;
    dvd=Idx.*Rd;
    count=1;
    
    %%% Current calculation loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while max(abs((Idx-Idxx)./Idx))>1e-10;
        count=count+1;
        if count>500, break, end
        Idxx=Idx;
        dvg=(Idx.*Rs+dvg)/2;
        dvd=(Idx.*Rd+dvd)/2;
        dvds=dvg+dvd; % total drop from source to drain
        
        Vdsi=Vds-dvds;
        Vgsi=Vgs-dvg;
        Vbsi=Vbs-dvg;
        
        Vsint=Vs+Idx.*(Rs0*1e-4/W)*dir;
        Vdint=Vd-Idx.*(Rd0*1e-4/W)*dir;
        Vgsraw=type*(Vg-Vsint);
        Vgdraw=type*(Vg-Vdint);
        
        
        % correct Vgsi and Vbsi
        % Vcorr is computed using external Vbs and Vgs but internal Vdsi
        % Qinv and Qinv_corr are computed with uncorrected Vgs, Vbs and
        % corrected Vgs, Vbs respectively.
        Vtpcorr=Vt0+gamma.*(sqrt(phib-Vbs)-sqrt(phib))-Vdsi.*delta;
        eVgpre = exp((Vgs-Vtpcorr)/(alpha*phit*1.5));
        FFpre = 1./(1+eVgpre);
        ab=2*(1-0.99*FFpre).*phit;
        Vcorr=(1+2.0*delta)*(ab./2.0).*(exp((-Vdsi)./(ab)));
        Vgscorr=Vgs+Vcorr-dvg;
        Vbscorr=Vbs+Vcorr-dvg;
        
        Vt0bs=Vt0+gamma.*(sqrt(phib-Vbscorr)-sqrt(phib));
        Vt0bs0=Vt0+gamma.*(sqrt(phib-Vbsi)-sqrt(phib));
        
        Vtp=Vt0bs-Vdsi.*delta-0.5*aphit;
        Vtp0=Vt0bs0-Vdsi.*delta-0.5*aphit;
        
        eVg=exp((Vgscorr-Vtp)/(aphit));
        FF=1./(1+eVg);
        eVg0=exp((Vgsi-Vtp0)/(aphit));
        FF0=1./(1+eVg0); 
        
        n=n0+abs(nd*Vdsi);
        nphit = n*phit;
        Qref=Cg*nphit;
        eta=(Vgscorr-(Vt0bs-Vdsi.*delta-FF*aphit))./(nphit);
        Qinv_corr = Qref.*log(1+exp(eta));
        eta0=(Vgsi-(Vt0bs0-Vdsi.*delta-FF*aphit))./(nphit); % compute eta0 factor from uncorrected intrinsic Vgs and internal Vds.
        %FF instead of FF0gives smoother C's!
        Qinv = Qref.*log(1+exp(eta0));
        
        vx0=vxo;
        Vdsats=vx0.*Leff./mu;
        Vdsat=Vdsats.*(1-FF)+ phit*FF;
        Fsat=(abs(Vdsi)./Vdsat)./((1+(abs(Vdsi)./Vdsat).^beta).^(1./beta));
        v=vx0.*Fsat;
        Idx = (W.*Qinv_corr.*v + 1*Idxx)/2;
    end
    %%% Current, positive into terminal y  %%%%%%%%%%%%%%%%%%%%%%%%
    Id(len_bias,1)=type*dir.*Idx; % in A
    Idlog(len_bias,1)=log10(Id(len_bias,1));
    Vdsi_out(len_bias,1) = Vdsi;
    
    
    % BEGIN CHARGE MODEL
    Vgt=Qinv./Cg;
    
    %Approximate solution for psis is weak inversion
    psis=phib+(1-gamma)/(1+gamma)*phit.*(1+log(log(1+exp((eta0)))));
    a=1+gamma./(2*sqrt(psis-(Vbsi))); % body factor
    Vgta=Vgt./a;   % Vdsat in strong inversion
    Vdsatq=sqrt(FF0.*(alpha*phit).^2+(Vgta).^2);  % Vdsat approx. to extend to weak inversion;
    % The multiplier of phit has strong effect on Cgd discontinuity at Vds=0.
    
    %Modified Fsat for calculation of charge partitioning (DD-NVSAT)
    betaq=beta;
    Fsatq=(abs(Vdsi)./Vdsatq)./((1+(abs(Vdsi)./Vdsatq).^betaq).^(1./betaq));
    x=1-Fsatq;
    
    
    % MVS 1.0.1 release uses Tsividis model for implementing DD-NVSAT
    % charges. In MVS 1.1.0, we use the formulation from L. Wei's 2012 TED
    % paper referenced in the documentation. Both models produce identical results. 
    
    % *** From Tsividis (DD-NVSAT model)
    % den=15*(1+x).^2;
    % qsc=Qinv*(6+12*x+8*x.^2+4*x.^3)./den;
    % qdc=Qinv*(4+8*x+12*x.^2+6*x.^3)./den;
    % qi=qsc+qdc;

    
    % *** From L. Wei (DD-NVSAT model)
    A=(1-x).^2./(12*(1-(1-x)/2));
    B=(5-2*(1-x))./(10-(1-x)/2);
    qsc=Qinv.*(0.5-(1-x)/6+A.*(1-B));
    qdc=Qinv.*(0.5-(1-x)/3+A.*B);
    qi=qsc+qdc;
  
    % Calculation of "ballistic" channel charge partitioning factors, qsb and qdb.
    % Default --> channel potential is assumed to increase linearly from the
    % virtual source point. That is V(x) = Vdsi*(x/Leff). 
    % Commented option --> channel potential profile is assumed to increase
    % parabolically from the virtual source point. That is, V(x) =
    % Vdsi*(x/Leff)^2.
   
    if (Vds < 1e-3)
        kq2=2.*qe/me*(zeta*Vdsi)/(vx0*vx0)*1e4;
        kq4=kq2.*kq2;
        % *** uncomment following two lines for linear potential profile
        qsb = Qinv*(0.5 - kq2/12.0 + kq4/32.0);
        qdb = Qinv*(0.5 - kq2/6.0 + 3*kq4/32.0);
        
        % *** uncomment following two lines for parabolic potential profile
        %qsb=Qinv*(0.5 - kq2/24.0 + kq4/80.0);
        %qdb=Qinv*(0.5 - kq2/8.0 + kq4/16.0);
        
    else
        kq=sqrt(2.*qe./me.*(zeta*Vdsi))./vxo.*1e2;  % 1e2 to convert cm/s to m/s. kq is unitless
        kq2=kq.^2;
        % *** uncomment following two lines for linear potential profile
        qsb=Qinv.*(4*(kq2+1).*sqrt(kq2+1)-(6*kq2+4))./(3*kq2.*kq2); 
        qdb=Qinv.*(2*(kq2-2).*sqrt(kq2+1)+4)./(3*kq2.*kq2); 
        
        % *** uncomment following two lines for parabolic potential profile
        %qsb=Qinv.*(asinh(sqrt(kq2))./sqrt(kq2)-(sqrt(kq2+1)-1)./kq2); 
        %qdb=Qinv.*((sqrt(kq2+1)-1)./kq2); 
        
    end
    
    % Flag for classic or ballistic charge partitioning:
    if (CTM_select == 1)   % classic DD-NVSAT
        qs=qsc;
        qd=qdc;
    else % ballistic blended with classic D/D
        Fsatq2=Fsatq.^2;
        qs=qsc.*(1-Fsatq2)+qsb.*Fsatq2;
        qd=qdc.*(1-Fsatq2)+qdb.*Fsatq2;
    end
    
    % Body charge based on approximate surface potential (psis) calculation.
    % With delta=0 using psis=phib in Qb gives continuous Cgs, Cgd, Cdd in SI,
    % while Cdd is smooth anyway.
    Qb(len_bias,1)=-type*W*Leff*(Cg*gamma*sqrt(psis-Vbsi)+(a-1)./a.*Qinv.*(1-qi));
    
    % DIBL effect on drain charge calculation.
    % Calculate dQinv at virtual source due to DIBL only.  Then:
    % Correct the qd factor to reflect this channel charge change due to Vds
    etai=(Vgsi-(Vt0bs0-FF*aphit))./(nphit); % Vt0bs0 and FF=FF0 causes least
    %discontinuity in Cgs and Cgd but produces a spike in Cdd at Vds=0 (in
    %weak inversion.  But bad in strong inversion)
    Qinvi = Qref.*(log(1+exp(etai)));
    dQinv=Qinv-Qinvi;
    dibl_corr=(1-FF0).*(1-Fsatq).*qi.*dQinv./Qinv;
    qd=qd-dibl_corr; %Potential problem area!
    
    % Inversion charge partitioning to terminals s and d accounting for
    % source drain reversal.
    Qinvs=type*Leff.*((1+dir).*qs+(1-dir).*qd)/2;
    Qinvd=type*Leff.*((1-dir).*qs+(1+dir).*qd)/2;
    
    % Overlap and outer fringe S and D to G charges
    % First calculate internal Vx and Vy
    
    Qxov=Cofs*(Vg-Vsint);
    Qyov=Cofd*(Vg-Vdint);
    
    % Inner fringing S and D to G charges; both screened by inversion at
    % that terminal (via FF function)
    Vt0x=Vt0+gamma*(sqrt(phib-type*(Vb-Vsint))-sqrt(phib));
    Vt0y=Vt0+gamma*(sqrt(phib-type*(Vb-Vdint))-sqrt(phib));
    Fs=1+exp((Vgsraw-(Vt0x-Vdsi.*delta.*Fsat)+aphit*0.5)./(1.1*nphit));
    Fd=1+exp((Vgdraw-(Vt0y-Vdsi.*delta.*Fsat)+aphit*0.5)./(1.1*nphit));
    
    FFx=Vgsraw-nphit.*log(Fs);
    FFy=Vgdraw-nphit.*log(Fd);
    
    Qxif=(type*Cif+CC*Vgsraw).*FFx;
    Qyif=(type*Cif+CC*Vgdraw).*FFy;
    
    % Total charge at internal terminals x and y.
    Qs(len_bias,1)=-W*(Qinvs+Qxov+Qxif);
    Qd(len_bias,1)=-W*(Qinvd+Qyov+Qyif);
    
    % Final charge balance
    Qg(len_bias,1)=-(Qs(len_bias,1)+Qd(len_bias,1)+Qb(len_bias,1));
    
    
    %end
end
