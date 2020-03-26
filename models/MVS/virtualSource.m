function [ dvds, dvg, Idx ] = virtualSource( params,V,dir)
%virtualSource Computes via an iterative solver technique the voltages for
%dvds and dvg used to accurately compute Ids of the device. This code is
%vectorized so that it can compute these voltages in parallel for n devices
%of the same type (either nmos or pmos).


Vds = V(1,:);
Vgs = V(2,:);
Vbs = V(3,:);
Vdsi = V(4,:);
Vs = V(5,:);
Vd = V(6,:);
Vg = V(7,:);

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
    n = n0 + nd * Vdsi;                           % Total subthreshold swing factor taking punchthrough into account [unit-less]
    nphit = n * phit;                             % Product of n and phit [used as one variable]
    aphit = alpha * phit;                         % Product of alpha and phit [used as one variable]
    
    Cofs = 0*( 0.345e-12/ etov ) * dLg/ 2.0 + Cof;  % s-terminal outer fringing cap [F/cm]
    Cofd = 0*( 0.345e-12/ etov ) * dLg/ 2.0 + Cof;  % d-terminal outer fringing cap [F/cm]
    Leff = Lgdr - dLg;                            % Effective channel length [cm]. After subtracting overlap lengths on s and d side

%      if (any(Vds < 0))
%             x = 0;
%             Vds = abs(Vds);
%         end

    Vt0bs=Vt0+gamma*(sqrt(abs(phib-Vbs))-sqrt(phib));
    
    % Denormalize access resistances and allocate them the "source" and
    % drain according to current flow
    Rs = 1e-4./ W * Rs0;                     % s-terminal resistance [ohms]
    Rd = Rs;  
    
    n=n0+nd*Vds;
    aphit = alpha*phit;
    nphit = n*phit;
    Qref=Cg*nphit;

%%% Initial values for current calculation %%%%%%%%%%%%%%%%%%%%
    FF=1./(1+exp((Vgs-(Vt0bs-Vds.*delta-0.5*aphit))/(aphit)));
    Qinv_corr = Qref.*log1p(exp((Vgs-(Vt0bs-Vds.*delta-FF*aphit))./(nphit)));
    Qinv = Qref.*log1p(exp((Vgs-Vt0bs)./(nphit)));
    Rt=Rs + Rd + (Lgdr-dLg)./(W.*Qinv*mu);

    vx0=vxo;
    Vdsats=W.*Qinv.*vx0.*Rt;
    Vdsat=Vdsats.*(1-exp(-Qinv./Qref))+phit*exp(-Qinv./Qref);
    Fsat=(1-exp(-2*Vds./Vdsat))./(1+exp(-2*Vds./Vdsat));
    Idx = W.*Fsat.*Qinv_corr.*vx0;
    Idxx=1e-15;
    dvg=Idx.*Rs;
    dvd=Idx.*Rd;
    count=1;
    dvds = 0;
    excessiveCountFlag = false;
   
    
    %%% Current calculation loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Occasionally takes steps that are too big!!%%
    %% Want to explore using a better root finder...%
     % counter of how many times I have needed to scale back the overshooting
    while max(abs((Idx-Idxx)./Idx))>1e-10
        count=count+1;
        if count>500, break, end
        Idxx=Idx;
        dvg=(Idx.*Rs+dvg)/2;
        dvd=(Idx.*Rd+dvd)/2;
        dvds=dvg+dvd; % total drop from source to drain
       
        
        if(any(dvds>Vds))
            divCount = 1;
            while(any(dvds>Vds))
                index = find(dvds > Vds);
                dvds(index) = 0.5*dvds(index);
                divCount = divCount + 1;
                if divCount > 500,break, end %prevent an infinite loop
            end        
        end
        
        
        Vdsi=Vds-dvds;
        Vgsi=Vgs-dvg;
        Vbsi=Vbs-dvg;
        
         if (any(Vdsi<0))
             x = 0;
        end
        
%         Vsint=Vs+Idx.*(Rs0*1e-4./W).*dir;
%         Vdint=Vd-Idx.*(Rd0*1e-4./W).*dir;
%         Vgsraw=type*(Vg-Vsint);
%         Vgdraw=type*(Vg-Vdint);
        
        
        % correct Vgsi and Vbsi
        % Vcorr is computed using external Vbs and Vgs but internal Vdsi
        % Qinv and Qinv_corr are computed with uncorrected Vgs, Vbs and
        % corrected Vgs, Vbs respectively.
        Vtpcorr=Vt0+gamma.*(sqrt(abs(phib-Vbs))-sqrt(phib))-Vdsi.*delta;
        eVgpre = exp((Vgs-Vtpcorr)/(alpha*phit*1.5));
        FFpre = 1./(1+eVgpre);
        ab=2*(1-0.99*FFpre).*phit;
        Vcorr=(1+2.0*delta)*(ab./2.0).*(exp((-Vdsi)./(ab)));
        Vgscorr=Vgs+Vcorr-dvg;
        Vbscorr=Vbs+Vcorr-dvg;
        
        Vt0bs=Vt0+gamma.*(sqrt(abs(phib-Vbscorr))-sqrt(phib));
        Vt0bs0=Vt0+gamma.*(sqrt(abs(phib-Vbsi))-sqrt(phib));
        
        
        
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
        Qinv_corr = Qref.*log1p(exp(eta));
        eta0=(Vgsi-(Vt0bs0-Vdsi.*delta-FF*aphit))./(nphit); % compute eta0 factor from uncorrected intrinsic Vgs and internal Vds.
        %FF instead of FF0gives smoother C's!
        Qinv = Qref.*log1p(exp(eta0));
        
        vx0=vxo;
        Vdsats=vx0.*Leff./mu;
        Vdsat=Vdsats.*(1-FF)+ phit*FF;
        Fsat=(Vdsi./Vdsat)./((1+(Vdsi./Vdsat).^beta).^(1./beta));
        v=vx0.*Fsat;
        Idx = ((W.*Qinv_corr.*v + 1*Idxx)/2);
        
        if(any(~isreal(Fsat)))
            warning('imaginary components in Fsat')
        end
        
        if( any(isnan(Idx)) || any(isinf(Idx)))
           error('Idx has NaN or Inf value')
        end
         
        if(count > 250)
            excessiveCountFlag = true;
        end
    end
    
    if(excessiveCountFlag)
     %   warning('count exceeded 450')
    end
    

end

