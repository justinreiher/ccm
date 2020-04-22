% This is a script which stages functions used to analyze the difference
% between the benchmark and offset designs. See
% offsetBenchTwoFlopComparePlots.m to get the results.

%voltage sources
gnd = vsrc('gnd',0,1);
vdd = vsrc('vdd',1,1);
clk = vtanhClock('clk',1e10,-pi/2,0.5,0.5,4,3e-10,0,1); 
clkbar = vtanhClock('clkbar',1e10,pi/2,0.5,0.5,4,3e-10,1,1);
din = vtanhDelay('din',0.5,0.5,5e10,true);

transistorWidths = ones(42,1)*450e-7;
syncAnalysisBench = TWO_FLOP_PG_SYNC('syncBenchAnalysis',transistorWidths,1);

nbAnalysisBench = nestedBisectionAnalysis(syncAnalysisBench,{syncAnalysisBench.vdd,syncAnalysisBench.gnd,syncAnalysisBench.clk,syncAnalysisBench.clkbar,syncAnalysisBench.d},...
   {vdd,gnd,clk,clkbar,vdd},{syncAnalysisBench.d,vdd},clk,7,'EKV','PTM 45nmHP');

transistorWidthsOffset = transistorWidths;
transistorWidthsOffset(8) = 600e-7;
transistorWidthsOffset(7) = 300e-7;

syncAnalysisOff = TWO_FLOP_PG_SYNC('synOffsetAnalysis',transistorWidthsOffset,1);
nbAnalysisOff = nestedBisectionAnalysis(syncAnalysisOff,{syncAnalysisOff.vdd,syncAnalysisOff.gnd,syncAnalysisOff.clk,syncAnalysisOff.clkbar,syncAnalysisOff.d},...
   {vdd,gnd,clk,clkbar,vdd},{syncAnalysisOff.d,vdd},clk,7,'EKV','PTM 45nmHP');

syncRef = TWO_FLOP_PG_SYNC('sync',1,1);
numDevNodes = 4*length(transistorWidths);
numStates = syncRef.nodeNum;


mask = [syncRef.cctPath.PGFF_0.PGL_0.q,...
        syncRef.cctPath.PGFF_0.q,...
        syncRef.cctPath.PGFF_1.PGL_0.q,...
        syncRef.cctPath.q];

uNodes = [syncRef.cctPath.PGFF_1.PGL_0.q,...
          syncRef.cctPath.PGFF_1.PGL_0.x0];

uVcrit = zeros(1,syncRef.nodeNum);
uVcrit(uNodes) = [1,-1];

span = [0 1.3e-9];
tCrit = 1.2e-9;


try
    load('benchAnalysisEKV')
    load('benchAnalysisEKV_RawDerivativesAll')
    dIdxDevAll = pchip(t,dIdx_Dev);
    CdevAll    = pchip(t,C_Dev);
    dIdxDev_t  = @(t) reshape(ppval(dIdxDevAll,t),[numDevNodes,numStates]);
    Cdev_t     = @(t) reshape(ppval(CdevAll,t),[numDevNodes,numStates]);
catch
    nbBench = nestedBisection(syncAnalysisBench,{syncAnalysisBench.vdd,syncAnalysisBench.gnd,syncAnalysisBench.clk,syncAnalysisBench.clkbar,syncAnalysisBench.d},...
    {vdd,gnd,clk,clkbar,vdd},[5.5e-10 4e-10],{syncAnalysisBench.d,din},clk,mask,'EKV','PTM 45nmHP');
    
    benchRuns = nbBench.simulate('2flopSyncBenchEKV',span,tCrit);
    [tbench,lambdaBench] = nbAnalysisBench.bisectionAnalysis('benchAnalysisEKV',benchRuns,uVcrit,tCrit);
    load('benchAnalysisEKV')
    load('benchAnalysisEKV_RawDerivativesAll')
    dIdxDevAll = pchip(t,dIdx_Dev);
    CdevAll    = pchip(t,C_Dev);
    dIdxDev_t  = @(t) reshape(ppval(dIdxDevAll,t),[numDevNodes,numStates]);
    Cdev_t     = @(t) reshape(ppval(CdevAll,t),[numDevNodes,numStates]);
end

g = zeros(1,length(t));
for i = 1:length(t)
    g(i) = wNorm_t(t(i))'*beta_t(t(i));
end
g = pchip(t,g);
g_t = @(t) ppval(g,t);
benchV = splMeta(t);
bench = pchip(t,benchV);
benchV = @(t) ppval(bench,t);
benchU = wNorm_t;

invGBench = @(in,out,span) componentGain(in,out,dIdxDev_t,Cdev_t,wNorm_t,g_t,span,nbAnalysisBench);
txGBench   = @(in,out,span,type) componentTxGain(in,out,dIdxDev_t,Cdev_t,wNorm_t,type,span,nbAnalysisBench);
gammaBench = @(in,out,t) gamma(in,out,t,map,dIdxDev_t,Cdev_t,wNorm_t,g_t,nbAnalysisBench);

JdevBench  = @(in,out,t) jacobianDev(in,out,t,dIdxDev_t,Cdev_t,nbAnalysis);

try
    load('offsetAnalysisEKV')
    load('offsetAnalysisEKV_RawDerivativesAll')
    dIdxDevAll = pchip(t,dIdx_Dev);
    CdevAll    = pchip(t,C_Dev);
    dIdxDev_t  = @(t) reshape(ppval(dIdxDevAll,t),[numDevNodes,numStates]);
    Cdev_t     = @(t) reshape(ppval(CdevAll,t),[numDevNodes,numStates]);
catch
    nbOffset = nestedBisection(syncAnalysisOff,{syncAnalysisOff.vdd,syncAnalysisOff.gnd,syncAnalysisOff.clk,syncAnalysisOff.clkbar,syncAnalysisOff.d},...
    {vdd,gnd,clk,clkbar,vdd},[5.5e-10 4e-10],{syncAnalysisOff.d,din},clk,mask,'EKV','PTM 45nmHP');
    offsetRuns = nbOffset.simulate('2flopSyncOffsetEKV',span,tCrit);
    [toffset,lambdaOffset] = nbAnalysisOff.bisectionAnalysis('offsetAnalysisEKV',offsetRuns,uVcrit,tCrit);
    load('offsetAnalysisEKV')
    load('offsetAnalysisEKV_RawDerivativesAll')
    dIdxDevAll = pchip(t,dIdx_Dev);
    CdevAll    = pchip(t,C_Dev);
    dIdxDev_t  = @(t) reshape(ppval(dIdxDevAll,t),[numDevNodes,numStates]);
    Cdev_t     = @(t) reshape(ppval(CdevAll,t),[numDevNodes,numStates]);
end


g = zeros(1,length(t));
for i = 1:length(t)
    g(i) = wNorm_t(t(i))'*beta_t(t(i));
end
g = pchip(t,g);
g_t = @(t) ppval(g,t);
offsetV = splMeta(t);
offset = pchip(t,offsetV);
offsetV = @(t) ppval(offset,t);
offsetU = wNorm_t;



invGOffset = @(in,out,span) componentGain(in,out,dIdxDev_t,Cdev_t,wNorm_t,g_t,span,nbAnalysisOff);
txGOffset   = @(in,out,span,type) componentTxGain(in,out,dIdxDev_t,Cdev_t,wNorm_t,type,span,nbAnalysisOff);
gammaOffset = @(in,out,t) gamma(in,out,t,dIdxDev_t,Cdev_t,wNorm_t,g_t,nbAnalysisOff);

nodePlotBench = @(in,out,t,g) plotGain(in,out,t,g,benchV,clk);
nodePlotOffset = @(in,out,t,g) plotGain(in,out,t,g,offsetV,clk);

JdevOffset = @(in,out,t) jacobianDev(in,out,t,dIdxDev_t,Cdev_t,nbAnalysisOff);

JdevPlot = @(in,out,start,stop) plotJacobianElementsAndU(in,out,start,stop,benchV,offsetV,benchU,offsetU,clk,JdevOffset,JdevBench);

componentCompare = @(in,out,tb,gb,to,go,legendStrings) compareGain(in,out,tb,gb,to,go,benchV,offsetV,benchU,offsetU,clk,legendStrings);

function plotJacobianElementsAndU(in,out,start,stop,benchV,offsetV,benchU,offsetU,clk,JdevOffset,JdevBench)
    
    t = linspace(start,stop,1000);
    vB = benchV(t);
    vO = offsetV(t);
    uB = benchU(t);
    uO = offsetU(t);
    
    jTemplateOff = JdevOffset(in,out,1);
    jTemplateOff = jTemplateOff(6:end,6:end);
    
    jTemplateBench = JdevBench(in,out,1);
    jTemplateBench = jTemplateBench(6:end,6:end);
    
    jDevTotOffset = zeros(length(reshape(jTemplateOff,[],1)),length(t));
    jDevTotBench  = zeros(length(reshape(jTemplateBench,[],1)),length(t));
    
    for i = 1:length(t)
        jOff = JdevOffset(in,out,t(i));
        jOff  = jOff(6:end,6:end);
        
        jBench = JdevBench(in,out,t(i));
        jBench = jBench(6:end,6:end);
        
        jDevTotOffset(:,i) = reshape(jOff,[],1);
        jDevTotBench(:,i)  = reshape(jBench,[],1);
    end
    
    
    yyaxis left
    hold on
    plot(t,vB(in-5,:),'-b',t,vO(in-5,:),'b--','linewidth',1.2)
    plot(t,vB(out-5,:),'-r',t,vO(out-5,:),'r--','linewidth',1.2)
    plot(t,clk.V(t),'-g','linewidth',1.2)
    plot(t,abs(uB(in-5,:)),'c-',t,abs(uO(in-5,:)),'c--','linewidth',1.2)
    plot(t,abs(uB(out-5,:)),'-k',t,abs(uO(out-5,:)),'k--','linewidth',1.2)
    xlabel('Time [s]')
    ylabel('Voltage [V]')
    yyaxis right
    plot(t,jDevTotBench,'-m',t,jDevTotOffset,'--m','linewidth',1.2)
    ylabel('J_{dev}[1/s]')
    legend('bench_{in}','offset_{in}','bench_{out}','offset_{out}','clk','bench |u_{in}|','offset |u_{in}|',...
        'bench |u_{out}|','offset |u_{out}|')

end

function compareGain(in,out,tb,gb,to,go,benchV,offsetV,benchU,offsetU,clk,legendStrings)
    
    %%%
    % legendStrings - input node name, output node name and component name
    tbV = linspace(tb(1),tb(end),1000);
    toV = linspace(to(1),to(end),1000);
    vB = benchV(tbV);
    vO = offsetV(toV);
    uB = benchU(tbV);
    uO = offsetU(toV);
	
	midClk = @(t) clk.V(t) - 0.5;
	tOff = fzero(midClk,3e-10);
	green = [0,0.5,0];
    blue  = [0,0.4770,0.7410];
    magenta = [1,0,1];
    f = figure;
    set(f,'defaultAxesColorOrder',[blue;green]);
    subplot(3,1,1)
    hold on
    title(strcat(legendStrings{3},' Bench vs Offset Compare'),'interpreter','latex')
    plot(tbV-tOff,vB(in-5,:),'-b',toV-tOff,vO(in-5,:),'b--','linewidth',1.5)
    plot(tbV-tOff,vB(out-5,:),'-r',toV-tOff,vO(out-5,:),'r--','linewidth',1.5)
    plot(tbV-tOff,clk.V(tbV),'-','color',green,'linewidth',1.2)
    %plot(tbV-tOff,abs(uB(in-5,:)),'c-',toV-tOff,abs(uO(in-5,:)),'c--','linewidth',1.5)
    %plot(tbV-tOff,abs(uB(out-5,:)),'-k',toV-tOff,abs(uO(out-5,:)),'k--','linewidth',1.5)
    ylabel('Signals [V]')% and |u(t)|')
    legend(strcat(legendStrings{1},' - Bench'),strcat(legendStrings{1},' - Offset'),...
        strcat(legendStrings{2},' - Bench'),strcat(legendStrings{2}, ' - Offset'), ...
        'clk','interpreter','latex')%,...
        %strcat(legendStrings{1},' $|u(t)|$ - Bench'),strcat(legendStrings{1},' $|u(t)|$ - Offset'),...
        %strcat(legendStrings{2},' $|u(t)|$ - Bench'),strcat(legendStrings{2},' $|u(t)|$ - Offset'),'interpreter','latex')
   
    set(gca,'fontsize',18)
    axis([tbV(1)-tOff,tbV(end)-tOff,0,1])
    set(f,'defaultAxesColorOrder',[magenta;green]);
    subplot(3,1,2)
    hold on
    yyaxis left
    plot(tb-tOff,log10(abs(gb)),'-m',to-tOff,log10(abs(go)),'--m','linewidth',1.5)
    ylabel('log_{10} g(t)')
    yyaxis right
    ylabel('clk [V]')% and |u(t)|')
    plot(tbV-tOff,clk.V(tbV),'-','color',green,'linewidth',1.2)
    %plot(tbV-tOff,abs(uB(in-5,:)),'c-',toV-tOff,abs(uO(in-5,:)),'c--','linewidth',1.5)
    %plot(tbV-tOff,abs(uB(out-5,:)),'-k',toV-tOff,abs(uO(out-5,:)),'k--','linewidth',1.5)
    legend(strcat(legendStrings{3},' $\log_{10} g(t)$ - Bench'), strcat(legendStrings{3}, ' $\log_{10} g(t)$ - Offset'), 'clk','interpreter','latex')%,...
     %   strcat(legendStrings{1},' $|u(t)|$ - Bench'),strcat(legendStrings{1},' $|u(t)|$ - Offset'),...
     %   strcat(legendStrings{2},' $|u(t)|$ - Bench'),strcat(legendStrings{2},' $|u(t)|$ - Offset'),...
     %'interpreter','latex')
    set(gca,'fontsize',18)
    axis([tbV(1)-tOff,tbV(end)-tOff,-inf,inf])
    
    gB = pchip(tb,gb);
    gBench = @(t) ppval(gB,t);
    
    gO = pchip(to,go);
    gOffset = @(t) ppval(gO,t);
    
    tEnd = min(tb(end),to(end));
    tStart = max(tb(1),to(1));
    tNew = linspace(tStart,tEnd,500);
    set(f,'defaultAxesColorOrder',[blue;green]);
    subplot(3,1,3)
    hold on
    yyaxis left
    plot(tNew-tOff, log10(gOffset(tNew)./gBench(tNew)),'-','linewidth',1.5)
    ylabel('$\log_{10} \left( \frac{g_{\mathit{offset}}(t)}{g_{\mathit{bench}}(t)}\right)$','interpreter','latex')
    yyaxis right
    xlabel('Time [s]')
    ylabel('clk [V]')% and |u(t)|')
    plot(tNew-tOff,clk.V(tNew),'-','color',green,'linewidth',1.2)
    %plot(tbV-tOff,abs(uB(in-5,:)),'c-',toV-tOff,abs(uO(in-5,:)),'c--','linewidth',1.5)
    %plot(tbV-tOff,abs(uB(out-5,:)),'-k',toV-tOff,abs(uO(out-5,:)),'k--','linewidth',1.5)
    legend(strcat(legendStrings{3}, ' $\log_{10} \left( \frac{ g_{\mathit{offset}}(t)}{g_{\mathit{bench}}(t)} \right)$'),'clk','interpreter','latex')%,...
        %strcat(legendStrings{1},' $|u(t)|$ - Bench'),strcat(legendStrings{1},' $|u(t)|$ - Offset'),...
        %strcat(legendStrings{2},' $|u(t)|$ - Bench'),strcat(legendStrings{2},' $|u(t)|$ - Offset'),...
        %'interpreter','latex')
    set(gca,'fontsize',18)
    axis([tbV(1)-tOff,tbV(end)-tOff,-inf,inf])
end
    

function plotGain(in,out,t,g,v,clk)
    vNodes = v(t);
    figure
    hold on
    yyaxis left
    plot(t,vNodes(in-5,:),'-b',t,vNodes(out-5,:),'-r')
    plot(t,clk.V(t),'-g')
    ylabel('Voltage [V]')
    xlabel('Times [s]')
    yyaxis right
    plot(t,log10(g),'-m')
    ylabel('log_{10}(gain)')
    legend('input','output','clk','gain')
end
    
function [gdev,t] = componentGain(in,out,dIdx_t,C_t,wNorm_t,g_t,span,testBench)

    map = testBench.computeMapMatrixDevice(in,out);
    M = [testBench.computeMapMatrix{:}];
        
    odeSettings = odeset('RelTol',1e-6,'AbsTol',1e-6);
    [t,gdev] = ode45(@(t,g) f(M*C_t(t),map*dIdx_t(t),wNorm_t(t),g),span,1,odeSettings);
    
    C = cell(1,100);
    tCom = linspace(span(1),span(2));
    for i = 1:100
        Ctemp = M*C_t(tCom(i));
        C{i} = Ctemp(6:end,6:end);
    end
    save(strcat(testBench.circuit.name,'_Caps.mat'),'C')
    
%    gdev = exp(g);
end

function [gdev,t] = componentTxGain(in,out,dIdx_t,C_t,wNorm_t,type,span,testBench)
    map = testBench.computeMapMatrixDeviceTx(in,out,type);
    M   = [testBench.computeMapMatrix{:}];
    odeSettings = odeset('RelTol',1e-6,'AbsTol',1e-6);
    [t,g] = ode45(@(t,g) f(M*C_t(t),map*dIdx_t(t),wNorm_t(t),g),span,1,odeSettings);
    gdev = exp(g);
end

function dgdt = f(C,dIdx,u,g)
    C = C(6:end,6:end);
    dIdx = dIdx(6:end,6:end);
    jac = C\dIdx;   
    dgdt = u'*jac*u*g;
end

function gammaOut = gamma(in,out,t,dIdx,C,u,g,testBench)

    map = testBench.computeMapMatrixDevice(in,out);

    gammaOut = zeros(1,length(t));
    
    for i = 1:length(t)
        c = C(t(i));
        c = map*c;
        dIdx = dIdx(t(i));
        j = c\(map*dIdx);
        j = j(6:end,6:end);
        U = u(t(i));
        gammaOut(i) = U'*j*U*g(t(i));
    end
end

function Jdev = jacobianDev(in,out,t,dIdx,CDev,testBench)
    mapDev = testBench.computeMapMatrixDevice(in,out);
    C = mapDev*CDev(t);
    Jdev = C\(mapDev*dIdx(t));
end
        
    