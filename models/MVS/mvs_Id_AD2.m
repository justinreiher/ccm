function [Id_ret,J] = mvs_Id_AD2 ( params, biasVoltages) %Vd,Vg,Vs,Vb )
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
    
    %      if (any(Vds < 0))
%             x = 0;
%             Vds = abs(Vds);
%         end

    Vt0bs=Vt0+gamma*(sqrt(abs(phib-Vbs))-sqrt(phib));
    
    % Denormalize access resistances and allocate them the "source" and
    % drain according to current flow
    Rs = 1e-4./ W * Rs0;                     % s-terminal resistance [ohms]
    Rd = Rs;  
    
    kB = 8.617e-5;                                % Boltzmann constant [eV/K]
    phit = kB*Tjun;
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
    Ids0 = W.*Fsat.*Qinv_corr.*vx0;
    
    [~,~, Ids] = IdsSolver(params, [Vds;Vgs;Vbs;Vdsi;Vs;Vd;Vg],Ids0);
      % IdsAD = Ids .*dir.*type;
       
       mvsAD = gradientinit([Vd;Vg;Vs;Vb;Ids]);
       
       IdsAD = computeIdsAD(params,mvsAD);
    
    Id_ret = Ids;
    JTot = IdsAD.dx;
    J = -JTot(:,end).^(-1)*JTot(:,1:4);
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
% dV=Ids;        
%        
%      Vsint=Vs+Ids.*(Rs0*1e-4./W).*dir;
%      Vdint=Vd-Ids.*(Rd0*1e-4./W).*dir;
%      Vgsraw=type*(Vg-Vsint);
%      Vgdraw=type*(Vg-Vdint);
%     
%     
% 
%     % Third try
%     intermediate_value = struct('Vb', biasVoltages(4), 'LARGE_VALUE', LARGE_VALUE, 'SMALL_VALUE', SMALL_VALUE, 'Qinv', Qinv, 'eta0', eta0, 'phit', phit, 'Vsi', Vsint, 'Vdi', Vdint, 'Vbsi', Vbsi, 'Vdsi', Vdsi, 'Vgsi', Vgsi, 'Vgsraw', Vgsraw, 'Vgdraw', Vgdraw, 'FF0', FF0, 'aphit', aphit, 'nphit', nphit, 'me', me, 'vx0', vx0, 'Leff', Leff, 'Vt0bs0', Vt0bs0, 'FF', FF, 'Qref', Qref, 'dir', dir, 'Cofs', Cofs, 'Cofd', Cofd, 'Fsat', Fsat);
end

