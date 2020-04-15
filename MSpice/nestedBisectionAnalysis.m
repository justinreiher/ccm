% testbench properties
%   circuit: the circuit to be tested
%   sources: voltage sources for each circuit inputs
%   map:    mapping between voltage sources and circuit inputs
%
% constructor
%   testbench(circuit,ports,sources);
%     circuit: an object of "circuit" class.
%     ports:    an cell of "node" objects
%     sources: an cell of "circuit" objects which must implement the V(t) function.
%          the size of ports and sources must be the same.
%          the voltage source is applied to the corresponding port
% methods
%   It provides functions for simulation
%   [t,v,dv] = simulate(tspan,v0,opts): this function simulate the circuit.
%     tspan:  the range of simulate time
%     v0:     the initial circuit state
%     opts:   opt.integrator is a function handle for ODE integrator.
%             it is also used in the integrator function.
%     t:      simulated time point
%     v:      voltage of circuit nodes, each row for a time
%     dv:     derivatives of node voltages, each row for a time
%
%   dv_ode = dV_ode(t,v_ode): construct ODE for ODE nodes
%     t:      used to evaluates voltage of sources
%     v_ode:  voltage of ODE nodes
%     dv_ode: time derivative of ODE voltages
%
%   It also provides function for verification (only if voltage sources support verification interface)
%   [c,err] = dV_ldi(region,states): LDI model for the whole system
%     region: an object of "shape" class
%     states: each element is passed to each voltage source
%     c,err:  dv \IN c'*[V;1]+/-err.
%
%   [c,err] = dV_ldi_var(region,var,states): LDI model for one node
%
%   str = name: the class name of the circuit
%
%   n  = dim:   the number of dimensions of the system dyanamics
%
%   ind = find_port_index: wrapper of circuit's "find_port_index" funcition.
%
%   n = inputNum: the number of input sources
%
%   is = inputs: return the voltage sources (same as sources)
%
%

classdef nestedBisectionAnalysis < handle
    properties (GetAccess='public', SetAccess='private');
        % a circuit to be analysed
        synchronizerCCT =[];
        % voltage sources
        sources={};
        %the input source that was used as a variable for bisection even
        %though bisection is done on voltage state the input source needs a
        %variable with which the delay can vary, a requirement to the method.
        inputSource = [];
        clockSource = [];
        inputSourceIndex = 0;
        outputNode = 0;
        % mapping between circuit ports and voltage sources
        map=[];
        mapMatrix = [];
        %the model dictionary that knows about all the models to compute dv for
        %the various device types know for the simulator.
        modelDictionary = [];
        %sets the default model and the AD version required for the analysis.
        defaultModel = [];
        defaultModelAD = [];
        
        %sets simulation options, i.e. model accuracy
        tbOptions = [];
        
    end
    properties (GetAccess = 'protected', SetAccess = 'private')
        synchronizerCCTMap;
    end
    
    methods
        % to create a nestedBisection testbench you need to provide:
        % - synchronizerCCT: a synchronizer circuit description
        % - ports: the external input and output ports to the synchronizer
        % - sources: the sources that are connected to the inputs of the synchronizer
        % - dinInterval: an interval where on one end a data transition causes the
        %     synchronizer to settle high and one where it will settle low,
        %     of the form [tLow,tHigh];
        % - dinSource: the data transition source, one which allows to change when the
        %     data transision occrus, i.e. vtanhDelay and the port to which
        %     it is connected in the form {synchronizer.d,dinSource}
        % - clockSource: the voltage source used to generate the clock (<-maybe there is a
        %     better way to access this)
        % - timeWindow: the desired time window to determine how many bisection restarts
        %     will be required in Seconds.
        % - nodeMask: a set of voltage nodes to observe at clock transition edges. This
        %     mask has a one to one correspondance with the clock edges,
        %     which means that when a clock edge occurs the algorithm knows
        %     what latch output to look at to determine the final outcome of
        %     the synchronizer. The edges are automatically computed and
        %     determined by the provided settings.
        % - modelDictionary: a model dictionary which contains all the known models to the
        %     testbench is required.
        % OPTIONAL parameters, but if not provided have default behaviours
        % - defaultModel: the default model to use to simulate the synchronizer. The
        %     default is to use the MVS model.
        % - defaultModelProcess: the default model process to use in the simulation. The default
        %     is to use the fit PTM 45nmHP model card.
        % - varargin: the following can be set by providing the appropriate (key,value)
        %     pairs or a struct which contains all the required parameters.
        %     The set of (key,value) pairs that are required in the overall struct are as follows:
        %         - 'capModel': defines the capacitor model to use for
        %             simulation, by default all capacitors are connected to
        %             'gnd'
        %         - 'vdd': the vdd voltage rail for the circuit, default for
        %                  the MVS 45nmHP model card is 1.0V
        %         - 'temp': the global temperature for the simulation in
        %                   Kelvin, the default is 298, i.e. 25 degrees C
        %         - 'numParallelCCTs': the number of parallel circuits to run
        %                              in parallel during the simulation, default is 10.
        %         - 'integratorOptions': a struct which is filled in by
        %                                odeset, a stopping function is set by default which halts
        %                                the simulation when all voltage nodes are within the
        %                                specified digitalSimOptions.
        %         - 'digitalSimOptions': a struct which contains:
        %                                  - 'minSimTime': the time which the
        %                                  simulation runs before intervening
        %                                  to allow initial transients to settle
        %                                  and must be greater than the time
        %                                  interval in which the data
        %                                  transition which causes the
        %                                  synchrnizer to settle high occurs,
        %                                  i.e. if dinInterval = [1,3], then
        %                                  'minSimTime' > 3
        %                                 - 'threshold': is the threshold
        %                                 where all voltage nodes in the
        %                                 synchronizer must settle, default
        %                                 is 10% of the rails.
        %                                 - 'stopGain': is the multiplicative
        %                                 factor in the stopping function,
        %                                 default is 30.
        %         - 'clkEdgeSettings': is a struct which contains:
        %                                  - 'low': which is a percentage of
        %                                  vdd on the low end where the half
        %                                  way of clos is found, default is
        %                                  48%
        %                                  - 'high': which is a percentage of
        %                                  vdd on the high end where the half
        %                                  way of clock is found, default is
        %                                  52%
        %                                  - 'clkEdgeSpread': is how much the
        %                                  high and low can be apart from
        %                                  each other in seconds
        %                                  - 'numberOfEdges': is the number
        %                                  of clock edges to find, default is
        %                                  3. This finds both falling and
        %                                  rising edges.
        %         - 'tCritSettings':
        
        function this = nestedBisectionAnalysis(synchronizerCCT,ports,sources,dinSource,clockSource,metaResolutionNode,...
                modelDictionary,defaultModel,defaultModelProcess,varargin)
            if(nargin<1), error('not enough parameters'); end
            if(nargin<2), ports = {}; end
            if(nargin<3), sources = {}; end
            if(nargin<4), error('the input data source is required in the form {synchronizer.din,dinSource}');end
            if(nargin<5), error('the clock source needs to be provided');end
            if(nargin<6), error('the synchronizer output port needs to be provided');end
            if(nargin<7), error('a model dictionary is required'); end
            if(nargin<9), defaultModel = 'MVS'; defaultModelProcess = 'PTM 45nmHP'; end
            if(nargin<10)
                integratorOptions = odeset('RelTol',1e-6,'AbsTol',1e-8,...
                    'InitialStep',0.5e-12);
                integratorOptionsBeta = odeset('RelTol',1e-6,'AbsTol',1e-6,...
                    'InitialStep',0.5e-12);%,'Events',@this.betaRenorm);
                
                betaSimOptions = struct('minSimTime',6e-10,'threshold',1e6);
                
                clkEdgeSettings = struct('low',0.48,'high',0.52,'clkEdgeSpread',20e-12,'numberOfEdges',1);
                
                tCritSettings   = struct('clkToQ',65e-12,'margin',2,'threshold',1e-2);
                
                digitalSimOptions = struct('minSimTime',6e-10,'threshold',0.1,...
                    'stopGain',30);
                
          this.tbOptions = struct('capModel','full','capScale',1e10,'vdd',1.0,'temp',298,'numParallelCCTs',10,...
              'stepOut',1,'polyFitDegree',3,'polyFitError',0.5,'isAD',true,'transSettle',1e-10,...
              'plotOptions',false,'integratorOptions',integratorOptions,'clockEdgeSettings',clkEdgeSettings,...
              'tCritSettings',tCritSettings,'digitalSimOptions',digitalSimOptions,...
                    'betaSimOptions',betaSimOptions,'integratorOptionsBeta',integratorOptionsBeta);
            end
            
            if(length(ports)~=length(sources))
                error('incorrect size of sources and ports');
            end
            
            if(isempty(this.tbOptions))
                this.packageTestbenchOptions(varargin);
            end
            
            this.modelDictionary = modelDictionary;
            this.defaultModel = modelDictionary.getModel(defaultModel,defaultModelProcess,this.tbOptions);
            this.defaultModelAD = modelDictionary.getModel(strcat(defaultModel,' AD'),defaultModelProcess,this.tbOptions);
            c = synchronizerCCT;
            c = c.flatten.vectorize; % jr: circuit must be flattened and vectorized
            %if the defaultModel class cannot compute vectorized quantities it
            %needs to evaluate them one at time. Not the responsibility of the
            %testbench.
            this.synchronizerCCT = c;
            %  this.circuit = circuit;
            this.sources = sources;
            this.map = synchronizerCCT.find_port_index(ports);
            
            this.inputSource = dinSource{2};
            this.modelDictionary = modelDictionary();
            this.synchronizerCCTMap = this.synchronizerCCT.getCircuitMap;
            this.mapMatrix = this.computeMapMatrix;
            if(any(this.map==0))
                error('can not find ports');
            end
            
            % Compute the stopping condition for the digital simulation
            vdd = this.tbOptions.vdd;
            this.clockSource = clockSource;
            this.outputNode = metaResolutionNode - length(this.sources);
            
            %find the index at which the input source is evaluated at
            this.inputSourceIndex = this.synchronizerCCT.find_port_index(dinSource{1});
            this.tbOptions.numNodes = this.synchronizerCCT.nodeNum;
            
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Functions for analysis
        
        % This function analyses the nested bisection data produced for a
        % nested bisection run. The user provides the critical time at which
        % the design needs to have settled by, the unit vector which selects
        % the nodes that are sensitive to the resolution of metastability.
        function [t,lambda] = bisectionAnalysis(this,filename,bisectionData,uVcrit,tCrit,opts)
            %get the last settle high and settle low data transitions and take
            %the mid point. (these should be the same because by this point in
            %the bisection algorithm when the data transition occured should be
            %numerically indistiguishable.
            if(ischar(filename)||isstring(filename))
                fprintf('starting bisection analysis and saving data to filename: %s ...\n',filename);
            else
                error('filename needs to be a string or array of characters')
            end
            
            if(nargin < 3), error('required to provide bisectionData'), end
            if(nargin < 4), error('required to provide vector to select critical criteria'), end
            if(length(uVcrit) == this.synchronizerCCT.nodeNum)
                iss = this.is_src;
                uVcrit = uVcrit(~iss);
                %normalize the unit vector if it is not already normalized.
                uVcrit = uVcrit/norm(uVcrit);
            end
            numStates = this.synchronizerCCT.nodeNum - length(this.sources);
            if(length(uVcrit) ~= numStates)
                error('incompatible uVcrit vector')
            end
            if(nargin < 5), error('need to specify t_crit'), end
            
            dataTransition = bisectionData{1,end};
            dataTransition = (dataTransition(1)+dataTransition(2))/2;
            
            time = bisectionData{10,2};
            t0 = time(1);
            beta0 = bisectionData{6,2}(:,1);
            [splMeta,t,teola] = this.splineBisection(bisectionData,tCrit,uVcrit);
            tspan = [t0,teola];
            
            
            assert(numStates == length(splMeta(1)));
            [~,numDevN] = size(this.mapMatrix{1});
            [totalNumStates,numDevP] = size(this.mapMatrix{2});
            numDev = numDevN+numDevP;
            
            
            % Integrator
            %set the default integrator options if none are specified
            if(nargin < 6)
                opts = this.tbOptions.integratorOptions;
                optsBeta = this.tbOptions.integratorOptionsBeta;
            end
            if(isfield('Integrator', opts))
                integrator = opts.Integrator;
            else
                integrator = @ode45;
            end
            
            indAnalysis = (t >= t0);
            t = t(indAnalysis);
            indAnalysis = (t <= teola);
            t = t(indAnalysis);
            
            iter = length(t);
            scale = this.tbOptions.capScale;
            jac = zeros([numStates^2,iter]);
            dIdx_Dev = zeros([numDev*totalNumStates,iter]);
            C_Dev = zeros([numDev*totalNumStates,iter]);
            vdot = zeros([numStates,iter]);
            dhda = zeros([numStates,iter]);
            fprintf('starting computation for Jac(t) and dh(t)/da...\n')
            tic()
            for i = 1:iter
                vSample = splMeta(t(i));
                [J,dhda(:,i),vdot(:,i),dIdxDev,CDev] = this.computeJacVdot(t(i),dataTransition,vSample);
                jac(:,i) = reshape(J,[numStates^2,1]);
                dIdx_Dev(:,i) = reshape(dIdxDev,[numDev*totalNumStates,1]);
                C_Dev(:,i) = reshape(CDev,[numDev*totalNumStates,1]);
            end
            toc()
            
           save(strcat(filename,'_jacAll.mat'),'jac','dhda','vdot','dIdx_Dev','C_Dev');
%             
%             %only use 1/10th the data
%             t = t(1:10:end);
%             jac = jac(:,1:10:end);
%             jacDev = jacDev(:,1:10:end);
%             dhda = dhda(:,1:10:end);
            
            jacP = pchip(t,jac);
            dIdxP = pchip(t,dIdx_Dev);
            CdevP   = pchip(t,C_Dev);
            
            
 %           dIdxDev_t = @(t) reshape(ppval(dIdxP,t),[numDev,totalNumStates]);
 %           Cdev_t    = @(t) reshape(ppval(CdevP,t),[numDev,totalNumStates]);
            jac_t = @(t) reshape(ppval(jacP,t),[numStates,numStates]);
            
            dhdaP = pchip(t,dhda);
            dhda_t = @(t) ppval(dhdaP,t);
            
            startTime = t0;
            stopTime = tspan(2);
            
            fprintf('starting computation for beta(t)...\n')
            tic()
            [tb, beta_ode] = integrator(@(tb, beta_ode)(this.dbeta_ode(jac_t(tb),dhda_t(tb),beta_ode)), [startTime stopTime], beta0, optsBeta);
            beta_ode = beta_ode';
            
            beta = spline(tb,beta_ode);
            beta_t = @(t) ppval(beta,t);
            toc()
            
            fprintf('starting computation for w(t) = ||uVcrit*S(t)||...\n');
            tic()
            [tw,w_ode] = integrator(@(tw,w_ode)(this.dwdt_ode(jac_t(tw),w_ode)), [stopTime startTime],uVcrit,opts);
            toc()
            w = w_ode';
            [~,iter] = size(w);
            wNorm = zeros(numStates,iter);
            
            for i = 1:iter
                wNorm(:,i) = w(:,i)/norm(w(:,i));
            end
            
            wNormt = spline(tw,wNorm);
            wNorm_t = @(t) ppval(wNormt,t);
            
            [~,iter] = size(t);
            lambda = zeros(1,iter);
            gamma = zeros(1,iter);
            fprintf('starting computation for g(t) and dg(t)/dt\n');
            tic()
            for i = 1:iter
                t_sample = t(i);
                [g,dgdt,lambda(i)] = this.gainSync(wNorm_t(t_sample),jac_t(t_sample),...
                    beta_t(t_sample),dhda_t(t_sample));
            end
            toc()
            save(strcat(filename,'.mat'),'t','dataTransition','splMeta',...
                'beta_t','jac_t','dhda_t','wNorm_t','lambda','g','dgdt');
        end % analysis
        
        % This function is a wrapper of circuit.dV by setting input values
        % as source.V(t) for the Matlab integrator.
        function dbeta = dbeta_ode(~,Jac,dhda,beta_ode)
            dbeta =  Jac*beta_ode + dhda;
            
        end
        
        function [dGdw_tCrit,gamma,t] = dGdw(this,traj_t,jac_t,dhda_t,dataTransition,beta_t,timeSpan,uVcrit)
            
            if(length(uVcrit) == this.synchronizerCCT.nodeNum)
                iss = this.is_src;
                uVcrit = uVcrit(~iss);
                %ensures normalization
                uVcrit = uVcrit/norm(uVcrit);
            end
            numStates = this.synchronizerCCT.nodeNum - length(this.sources);
            if(length(uVcrit) ~= numStates)
                error('incompatible uVcrit vector')
            end

            deviceMap = this.synchronizerCCT.getCircuitMap;
            numUniqueElements = length(deviceMap);
            wid = cell(1,numUniqueElements);
            for i = 1:numUniqueElements
                wid{i} = ones(1,length(deviceMap{i}{2}));
            end
            
            [~,numAD_var] = populateADoptimizationVar(1,wid{1},wid{2});
            
            gamma0 = zeros(numStates,numAD_var);
            fprintf('Starting computation of dG/dw...\n') 
            tic()
            [t,gamma] = ode45(@(t,gamma) this.gamma_ode(t,traj_t(t),jac_t(t),dhda_t(t),beta_t(t),gamma,dataTransition),timeSpan,gamma0);
            dGdw_tCrit = expandW((uVcrit*reshape(gamma(end,:),numStates,numAD_var))');
            toc()
        end
        
        function gammaDot = gamma_ode(this,t,traj,jac,dhda,beta,gamma,dataTransition)
            
            vfull = this.vfull(t,dataTransition,traj);
            synchronizerCCTMapping = this.synchronizerCCTMap;
            [~,numUniqueElements] = size(synchronizerCCTMapping);

            
            % The following is all done to setup the AD variable
            Idev_f = cell(1,numUniqueElements);
            Ctot   = cell(1,numUniqueElements);
            Ms = this.mapMatrix;
            vfull(this.inputSourceIndex,:) = dataTransition;
            wid = cell(1,numUniqueElements);
            for i = 1:numUniqueElements
                %get the circuit element grouping
                elementGroup = synchronizerCCTMapping{i};
                %extract the device properties and associated map
                deviceProperties = elementGroup{1};
                deviceMap = elementGroup{2};
                %extract the device type and model along with any specific model
                %parameters.
                %   deviceModelParams = deviceProperties.deviceModelParams;
                %struct which contains relevant parameters for the device in a
                %one to one correspondance. I.e.:
                %the parameter with field w {width_device1, width_device2, ...}
                vectorizedParams = deviceProperties.parameters;
                %determine the number of terminals per device and the number of
                %devices. Figure out the number of points to evaluate and shape V
                %such that each column represents a device for a particular
                %point.
                [numDeviceTerminals,numDevices] = size(deviceMap);
                [numCircuitNodes,numPointsToEval] = size(vfull);
                
                wid{i} = [vectorizedParams.wid{:}];
            end
            widN = wid{1};
            widP = wid{2};

            %AD_f = hessianinit([vfull;widN';widP']);
            % The above allows all transistor widths to beADvariables, the
            % function populateADoptimizationVar(vfull,widN,widP) allows
            % the user to choose which transistors are allowed to be
            % uncoupled. 
            [AD_f,numAD_w] = populateADoptimizationVar(vfull,widN,widP);
            widPrev = [];
            
            for i = 1:numUniqueElements
                
                %get the circuit element grouping
                elementGroup = synchronizerCCTMapping{i};
                %extract the device properties and associated map
                deviceProperties = elementGroup{1};
                deviceMap = elementGroup{2};
                %extract the device type and model along with any specific model
                %parameters.
                deviceType = deviceProperties.deviceType;
                deviceModel = deviceProperties.deviceModel;
                %   deviceModelParams = deviceProperties.deviceModelParams;
                %struct which contains relevant parameters for the device in a
                %one to one correspondance. I.e.:
                %the parameter with field w {width_device1, width_device2, ...}
                vectorizedParams = deviceProperties.parameters;
                vectorizedParams.deviceMap = deviceMap;
                vectorizedParamsAD = vectorizedParams;
                %determine the number of terminals per device and the number of
                %devices. Figure out the number of points to evaluate and shape V
                %such that each column represents a device for a particular
                %point.
                wid = [vectorizedParams.wid{:}];
                [numDeviceTerminals,numDevices] = size(deviceMap);
                [numCircuitNodes,numPointsToEval] = size(vfull);
                 
                Vf_AD = AD_f(1:numCircuitNodes);
                widAD = cell(1,length(wid));
                for j = 1:length(wid)
                    widAD{j} = AD_f(numCircuitNodes+j+length(widPrev));
                end
                vectorizedParamsAD.wid = widAD;
                widPrev = wid;
                
                %set the voltage of the data transition node based on the time
                %and previously determined data transition time
                
                din = this.inputSource.V(t,Vf_AD(this.inputSourceIndex,:));
                Vf_AD(this.inputSourceIndex,:) = din;
                
                V = Vf_AD(deviceMap,:);
                
                % V(which_terminal,which_point_in_phase_space,which_device)
                V = reshape(V,[numDeviceTerminals,numPointsToEval,numDevices]);
                
                if(strcmp(deviceModel,'default'))
                    [Idev,capParams] = this.defaultModelAD.I(deviceType,vectorizedParamsAD,V,this.tbOptions);
                    vectorizedParamsAD.capParams = capParams;
                    CDevice = this.defaultModelAD.C(deviceType,vectorizedParamsAD,V,this.tbOptions);
                else
                    model = this.modelDictionary.getModel(deviceModel,deviceModelParams);
                    [Idev,capParams] = model.I(deviceType,vectorizedParamsAD,V,this.tbOptions);
                    vectorizedParams.capParams = capParams;
                    CDevice = this.model.C(deviceType,vectorizedParams,V,this.tbOptions);
                end
                
                Idev_f{i} = Ms{i}*Idev;
                Ctot{i} = Ms{i}*CDevice;
                
            end
        Idev = zeros(size(Idev_f{1}))*V(1);
        C    = zeros(size(Ctot{1}))*V(1);

        for i = 1:length(Idev_f)
            Idev = Idev + Idev_f{i};
            C    = C    + Ctot{i};
        end

        %compute currents for voltage sources if any exist and their model
        %follows the default
        for i = 1:numUniqueElements
            %get the circuit element grouping
            elementGroup = synchronizerCCTMapping{i};
            %extract the device properties and associated map
            deviceProperties = elementGroup{1};
            deviceMap = elementGroup{2};
            %extract the device type and model along with any specific model
            %parameters.
            deviceType = deviceProperties.deviceType;
            deviceModel = deviceProperties.deviceModel;
            %make the currents and capacitance on each side of the voltage
            %source the same.
            if(strcmp(deviceType,'vsrc') && strcmp(deviceModel,'default'))
                Idev(deviceMap,:) = Idev(deviceMap,:) + Idev(flipud(deviceMap),:);
            end
        end
        

        
        if(strcmp(this.tbOptions.capModel,'gnd'))
            
            isSource = this.is_src;
            I = Idev(~isSource);
            
            C = C(~isSource);
            f = I./C;
            
            offset = length(this.sources);
            totalNumNodes = this.synchronizerCCT.nodeNum;
            numStates = totalNumNodes - offset;
            alphaSource = this.inputSourceIndex;
            
            gamma = reshape(gamma,numStates,numAD_w);
            
            df_dwda = zeros(numStates,numAD_w);
            betaJac_fw = zeros(numStates,numAD_w);
            for i = 1:numStates
                df_dwda_i = f(i).hx;
                df_dwda(i,:) = df_dwda_i(totalNumNodes+1:end,alphaSource);
                betaJac_fw(i,:) = df_dwda_i(totalNumNodes+1:end,offset+1:totalNumNodes)*beta;
            end
            nonHomogeneous = df_dwda+betaJac_fw;
        else
            
            numDevices = length(widN) + length(widP);
            offset = length(this.sources);
            totalNumNodes = this.synchronizerCCT.nodeNum;
            numStates = totalNumNodes - offset;
            alphaSource = this.inputSourceIndex;
            
            isSource = this.is_src;
            I = Idev(~isSource);
            dI = I.dx;
            dIdx = dI(:,1+offset:totalNumNodes);
            dIdw = dI(:,1+totalNumNodes:end);
            dIda = dI(:,alphaSource);
            %C here is the diagonal matrices of the sparse matrix I want to have,
            %MATLAB doesn't offer a convienient way to create a block
            %diagonal sparse matrix other than by the following:
            B = mat2cell(sparse(C),numCircuitNodes,numCircuitNodes);
            C = blkdiag(B{:});
            C = C(~isSource,~isSource);
            Cap = C.x;
            dC  = C(:).dx; %all the first derivatives for C
            dCda = reshape(dC(:,alphaSource),numStates,numStates);
            f = Cap\I.x;
            Cinv_dIdx_Beta = Cap\(dIdx*beta);
            Cinv_dIdw      = Cap\dIdw;
            Cinv_dIda      = Cap\dIda;
            dCda_Cinv_dIdw = dCda*Cinv_dIdw;
            Cinv_dCda_f    = Cap\(dCda*f);
          
            
            %let try to be clever about this later, ouch though triple
            %loop...
            
            dCdwdx_f = zeros(numStates,numStates,numAD_w);
            dCdwda = zeros(numStates,numStates,numAD_w);
            
            for i = 1:numStates
                for j = 1:numStates
                    capElement = C(i,j).hx;
                    dCdwdx_f(i,j,:) = capElement(1+totalNumNodes:end,1+offset:totalNumNodes)*f;
                    dCdwda(i,j,:) = capElement(1+totalNumNodes:end,alphaSource);
                end
            end
            
            dCdwdx_f_Beta = zeros(numStates,numAD_w);
            dCdwda_f = zeros(numStates,numAD_w);

            for i = 1:numAD_w
                dCdwdx_f_Beta(:,i) = dCdwdx_f(:,:,i)*beta;
                dCdwda_f(:,i)      = dCdwda(:,:,i)*f;
            end
            %%% I think the above is correct
            
            gamma = reshape(gamma,numStates,numAD_w);
            
            dCdx_f = zeros(numStates,numStates);
            dIdwdx_Beta = zeros(numStates,numAD_w);
            dIdwda      = zeros(numStates,numAD_w);
            for i = 1:numStates
                dC_i = reshape(dC(:,i+offset),numStates,numStates);
                dCdx_f(i,:) = dC_i*f;
                ddI_i = I(i).hx;
                dIdwdx_Beta(i,:) = ddI_i(1+totalNumNodes:end,1+offset:totalNumNodes)*beta;
                dIdwda(i,:) = ddI_i(1+totalNumNodes:end,alphaSource);
            end
            
            Cinv_dCdx_f_Beta = Cap\(dCdx_f*beta);
            
            dCdw_f = zeros(numStates,numAD_w);
            dCdw_Cinv_dCdx_f_Beta = zeros(numStates,numAD_w);
            dCdw_Cinv_dIdx_Beta  = zeros(numStates,numAD_w);
            dCdw_Cinv_dCda_f     = zeros(numStates,numAD_w);
            dCdw_Cinv_dIda       = zeros(numStates,numAD_w);
            for i = 1:numAD_w
                dC_i = reshape(dC(:,i+totalNumNodes),numStates,numStates);
                dCdw_f(:,i) = dC_i*f;
                dCdw_Cinv_dCdx_f_Beta(:,i) = dC_i*Cinv_dCdx_f_Beta;
                dCdw_Cinv_dIdx_Beta(:,i)   = dC_i*Cinv_dIdx_Beta;
                dCdw_Cinv_dCda_f(:,i)      = dC_i*Cinv_dCda_f;
                dCdw_Cinv_dIda(:,i)        = dC_i*Cinv_dIda;
            end
            
            Cinv_dCdw_f = Cap \dCdw_f;
            dCda_Cinv_dCdw_f = dCda*Cinv_dCdw_f;
            dCdx_Cinv_dCdw_f_Beta = zeros(numStates,numAD_w);
            dCdx_Cinv_dIdw_Beta   = zeros(numStates,numAD_w);
            
            for i = 1:numStates
                dC_i = reshape(dC(:,i+offset),numStates,numStates);
                dCdx_Cinv_dCdw_f_Beta(i,:) = beta'*dC_i*Cinv_dCdw_f;
                dCdx_Cinv_dIdw_Beta(i,:)   = beta'*dC_i*Cinv_dIdw;
            end
            
            dJdw_Beta = (dCdw_Cinv_dCdx_f_Beta - dCdwdx_f_Beta + dCdx_Cinv_dCdw_f_Beta - dCdx_Cinv_dIdw_Beta - dCdw_Cinv_dIdx_Beta + dIdwdx_Beta);
            dfdwda    = (dCdw_Cinv_dCda_f - dCdwda_f + dCda_Cinv_dCdw_f - dCda_Cinv_dIdw - dCdw_Cinv_dIda + dIdwda);
            
            
            nonHomogeneous = Cap\(dJdw_Beta + dfdwda);
            
            
%             
%             df_dwda = zeros(numStates,numDevices);
%             betaJac_fw = zeros(numStates,numDevices);
%             dCdxf = zeros(numStates,numStates);
% 
%             for i = 1:numStates
%                 
%                 row_Ihx = I(i).hx;
%                 d2Idwdx = row_Ihx(offset+1:totalNumNodes,totalNumNodes+1:end);
%                 dC_i = dC(:,i+offset);
%                 dCdw_i = dC_i(:,totalNumNodes+1:end);
%                 dCdwJ_i = (dCdw_i'*jac)';
%                 betaJac_fw(i,:) = beta'*(-dCdwJ_i + d2Idwdx);
%                 
%                 row_Ihx = I(i).hx;
%                 
%                 dCdwDhda_i = -dhda'*dCdw_i;
%                 d2Idwda_i = row_Ihx(totalNumNodes+1:end,alphaSource)';
%                 
%                 df_dwda(i,:) = dCdwDhda_i + d2Idwda_i;
%             end
% 
%             nonHomogeneous = Cap\(betaJac_fw + df_dwda);
        end
        
        gammaDot = jac*gamma + nonHomogeneous;
        gammaDot = reshape(gammaDot,[],1);
        
    end
    
        function packageTestbenchOptions(this,tbOptions)
            if(isstruct(tbOptions{1}))
                this.tbOptions = tbOptions{1};
            else
                paramNum = length(tbOptions);
                if(~(mod(paramNum,2) == 0))
                    error('Simulation Options have (key,value) pair mismatch');
                end
                for i = 1:(paramNum/2)
                    index = 2*i-1;
                    key = tbOptions{index};
                    value = tbOptions{index+1};
                    options.(key) = value;
                end
                this.tbOptions = options;
            end
            
        end
        
        
        function [VmSpl_t,time,teola] = splineBisection(this,bisectionData,tCrit,uVcrit)
            %this function takes in a set of bisection trajectories and creates
            %a spline fit for the trajectory that is metastable up to time
            %teola
            
            % first we find trajectories which went just a bit too high and
            % too low of our settling criteria to find at tCrit a
            % trajectory which just fails.
            numParallelCCTs = this.tbOptions.numParallelCCTs;
            numNodesPerCCT  = this.tbOptions.numNodes - length(this.sources);
            output = this.outputNode;
            alphaRef = linspace(0,1,numParallelCCTs);
            deltaV = this.tbOptions.tCritSettings.threshold;
            
            
            %%%%%%
            % The following is specifically implemented for a simple
            % stopping criteria - FUTURE WORK is to make this more generic
            
            digSimSettleCriteria = this.tbOptions.digitalSimOptions.threshold;
            vdd = this.tbOptions.vdd;
            
            failLow = digSimSettleCriteria*vdd;
            failHigh = (1-digSimSettleCriteria)*vdd;
            
            VallEnd = bisectionData{14,end};
            tEnd  = bisectionData{10,end};
            
            ind = find(tEnd >= tCrit);
            ind = ind(1);
            
            vHind = 0;
            vLind = 1;
            
            for i = 1:this.tbOptions.numParallelCCTs
                ViQ = VallEnd(ind,output + (i-1)*numNodesPerCCT);
                ViQEnd = VallEnd(end,output + (i-1)*numNodesPerCCT);
                if( (ViQ < failHigh) && (ViQEnd > failHigh))
                    indH = i;
                    VhCritH = VallEnd(1,1+(indH-1)*numNodesPerCCT:indH*numNodesPerCCT)';
                end
                if( (ViQ > failLow) && (ViQEnd < failLow))
                    indL = i;
                    VlCritL = VallEnd(1,1+(indL-1)*numNodesPerCCT:indL*numNodesPerCCT)';
                end
            end
            
            figure
            hold on
            
            VH = VallEnd(:,1+(indH-1)*numNodesPerCCT:indH*numNodesPerCCT);
            VL = VallEnd(:,1+(indL-1)*numNodesPerCCT:indL*numNodesPerCCT);
            
            plot(tEnd,VH(:,output),'r')
            plot(tEnd,VL(:,output),'b')
            plot(tEnd(1),VhCritH(output),'rx')
            plot(tEnd(1),VlCritL(output),'bx')
            plot(tCrit*ones(1,100),linspace(0,1),'--')
            
            vh = pchip(tEnd,VH(:,output));
            vl = pchip(tEnd,VL(:,output));
            
            vhCrit = @(t) ppval(vh,t) - failHigh;
            vlCrit = @(t) ppval(vl,t) - failLow;
            
            indH = find( VH(:,output) - failHigh > 0);
            indL = find( VL(:,output) - failLow < 0 );
            
            tH = tEnd(indH(1));
            tL = tEnd(indL(1));
            
            % tH = fzero(vhCrit,tEnd(end));
            % tL = fzero(vlCrit,tEnd(end));
            
            plot(tH*ones(1,100),linspace(0,1),'r--')
            plot(tL*ones(1,100),linspace(0,1),'b--')
            
            delta_tH = tH - tCrit;
            delta_tL = tL - tCrit;
            
            
            indexPrev = bisectionData{7,end-1};
            VallPrev  = bisectionData{14,end-1};
            tPrev     = bisectionData{10,end-1};
            
            vEndPrev = bisectionData{4,end-2};
            vHEndPrev = vEndPrev(:,1);
            vLEndPrev = vEndPrev(:,2);
            
            VHprev = VallPrev(:,1+(indexPrev(2)-1)*numNodesPerCCT:indexPrev(2)*numNodesPerCCT);
            VLprev = VallPrev(:,1+(indexPrev(1)-1)*numNodesPerCCT:indexPrev(1)*numNodesPerCCT);
            
            plot(tPrev,VHprev(:,output),'-.r')
            plot(tPrev,VLprev(:,output),'-.b')
            
            Vh = pchip(tPrev,VHprev');
            Vl = pchip(tPrev,VLprev');
            
            vH_interpNew = ppval(Vh,tEnd(1) - delta_tH);
            vL_interpNew = ppval(Vl,tEnd(1) - delta_tH);
            
%              ind = abs(vH_interpNew - vL_interpNew) > 2.5e-3;
              vH_interpNew = vH_interpNew(output);
              vL_interpNew = vL_interpNew(output);
              VhCritH      = VhCritH(output);
            
            plot(tEnd(1) - delta_tH,vH_interpNew,'rs',tEnd(1) - delta_tH,vL_interpNew,'bs')
            
            alphaH = (vH_interpNew - vL_interpNew)'*(VhCritH-vL_interpNew)/((vH_interpNew - vL_interpNew)'*(vH_interpNew - vL_interpNew));
            pH =  (alphaH*vH_interpNew + (1-alphaH)*vL_interpNew);
            
            plot(tEnd(1)-delta_tH,pH,'ro')
            
            vH_interpNew = ppval(Vh,tEnd(1) - delta_tL);
            vL_interpNew = ppval(Vl,tEnd(1) - delta_tL);
            
%              ind = abs(vH_interpNew - vL_interpNew) < 2.5e-3;
              vH_interpNew = vH_interpNew(output);
              vL_interpNew = vL_interpNew(output);
              VlCritL      = VlCritL(output);
            
            alphaL = (vH_interpNew - vL_interpNew)'*(VlCritL-vL_interpNew)/((vH_interpNew - vL_interpNew)'*(vH_interpNew - vL_interpNew));
            pL = (alphaL*vH_interpNew + (1-alphaL)*vL_interpNew);
            
            plot(tEnd(1)-delta_tL,pL,'bo')
            
            alphaH = alphaH*alphaRef(indexPrev(2))+(1-alphaH)*alphaRef(indexPrev(1));
            alphaL = alphaL*alphaRef(indexPrev(2))+(1-alphaL)*alphaRef(indexPrev(1));
            
            newVHStart = alphaH*vHEndPrev + (1-alphaH)*vLEndPrev;
            newVLStart = alphaL*vHEndPrev + (1-alphaL)*vLEndPrev;
            
            metaStart = (alphaH + alphaL)/2*vHEndPrev + (1-(alphaH + alphaL)/2)*vLEndPrev;
            
            plot(tPrev(1),newVHStart(output),'ro')
            plot(tPrev(1),newVLStart(output),'bo')
            plot(tPrev(1),vHEndPrev(output),'r^')
            plot(tPrev(1),vLEndPrev(output),'b^')
            plot(tPrev(1),metaStart(output),'xm')
            
            teola = tPrev(1);

%             for i = 1:numParallelCCTs-1
%                 ViQ = VallEnd(ind,output+(i-1)*numNodesPerCCT);
%                 ViQEnd = VallEnd(end,output+ (i-1)*numNodesPerCCT);
%                 
%                 if( (ViQ < failHigh) && (ViQEnd > failHigh))
%                     if(ViQ > vHind)
%                         vHind = ViQ;
%                         for j = 1:numParallelCCTs
%                             Vh = VallEnd(ind,output+(j-1)*numNodesPerCCT);
%                             if(Vh > failHigh)
%                                 indH = [i,j];
%                                 break;
%                             end
%                         end
%                     end
%                 end
%                 
%                 if( (ViQ > failLow) && (ViQEnd < failLow))
%                     if(ViQ <vLind)
%                         vLind = ViQ;
%                         for j = 1:numParallelCCTs
%                             Vl = VallEnd(ind,output+(j-1)*numNodesPerCCT);
%                             if(Vl < failLow)
%                                 indL = [j,i];
%                                 break;
%                             end
%                         end
%                     end
%                 end
%             end
            
            %Now that we have those trajectories, we now linearly
            %interpolate between the set of trajectories to create a
            %trajectory which just fails to go high and one which just
            %fails to go low
            
%             VhH = pchip(tEnd',VallEnd(:,1+(indH(2)-1)*numNodesPerCCT:indH(2)*numNodesPerCCT)');
%             VhL = pchip(tEnd',VallEnd(:,1+(indH(1)-1)*numNodesPerCCT:indH(1)*numNodesPerCCT)');
%             
%             VlH = pchip(tEnd',VallEnd(:,1+(indL(2)-1)*numNodesPerCCT:indL(2)*numNodesPerCCT)');
%             VlL = pchip(tEnd',VallEnd(:,1+(indL(1)-1)*numNodesPerCCT:indL(1)*numNodesPerCCT)');
%             
%             VhHCrit = ppval(VhH,tCrit);
%             VhLCrit = ppval(VhL,tCrit);
%             
%             VlHCrit = ppval(VlH,tCrit);
%             VlLCrit = ppval(VlL,tCrit);
%             
%             alphaH = (failHigh - VhLCrit(output))./(VhHCrit(output)-VhLCrit(output));
%             alphaL = (failLow - VlLCrit(output))./(VlHCrit(output) - VlLCrit(output));
%             
%             VH = alphaH*VallEnd(:,1+(indH(2)-1)*numNodesPerCCT:indH(2)*numNodesPerCCT)+...
%                 (1-alphaH)*VallEnd(:,1+(indH(1)-1)*numNodesPerCCT:indH(1)*numNodesPerCCT);
%             
%             VL = alphaL*VallEnd(:,1+(indL(2)-1)*numNodesPerCCT:indL(2)*numNodesPerCCT)+...
%                 (1-alphaL)*VallEnd(:,1+(indL(1)-1)*numNodesPerCCT:indL(1)*numNodesPerCCT);
            
            %With the just failing high and just failing low trajectories,
            %we now compute teola by finding the larger time
            %uVcrit*VH - uVcrit*VL <= deltaV
            
            
%             uVHteola = @(t) ppval(uVH,t) - deltaV;
%             uVLteola = @(t) ppval(uVL,t) + deltaV;
%             
%             
%             
%             teolaH = fzero(uVHteola,tEnd(start));
%             teolaL = fzero(uVLteola,tEnd(start));
%             
%             teola = max(teolaH,teolaL);
            
%             uV = @(t) ppval(uVH,t) - ppval(uVL,t) - deltaV;
%             
%             teola = fzero(uV,tEnd(start));
            
            %With teola established, lets now spline the perfectly
            %metastable trajectory starting with the before last nested
            %bisection trajectories.
            
            %find times that are greater than teola to reject them
            indTend = find(tEnd >= teola);
            
            highTraj = bisectionData{8,end-2};
            lowTraj = bisectionData{9,end-2};
            tEnd2   = bisectionData{10,end-2};
            
            VHlast = alphaH*highTraj + (1-alphaH)*lowTraj;
            VLlast = alphaL*highTraj + (1-alphaL)*lowTraj;
            
            %spline the trajectory which just failed high and low
            
            VHchip = pchip(tEnd,VH');
            VLchip = pchip(tEnd,VL');
            
            %initialize data to be saved for splining at the end
            
            alphaM = zeros(1,length(bisectionData)-2);
            Vmeta = cell(1,length(bisectionData)-2);
            time = cell(1,length(bisectionData)-2);

            
            aH = zeros(1,length(bisectionData));
            aL = zeros(1,length(bisectionData));
            vH = cell(1,length(bisectionData));
            vL = cell(1,length(bisectionData));
            timeFull = cell(1,length(bisectionData));
            
            %spline the before last trajectory to compute the starting
            %interpolatant alphaH and alphaL
            
%             VHlastInterp = pchip(tEnd2,VHlast);
%             VLlastInterp = pchip(tEnd2,VLlast);
%             
%             alphaH = norm(ppval(VHchip,teola) - ppval(VLlastInterp,teola))/norm(ppval(VHlastInterp,teola)-ppval(VLlastInterp,teola));
%             alphaL = norm(ppval(VLchip,teola) - ppval(VLlastInterp,teola))/norm(ppval(VHlastInterp,teola)-ppval(VLlastInterp,teola));
%             
            %start with the "perfectly metastable trajectory" in the middle
            %between alphaH and alphaL
            alphaM(end) = 0.5*(alphaH+alphaL);
            %Vmeta{end} = (1-alphaM(end))*VLlast + alphaM(end)*VHlast;
            
            % for diagnostic purposes, also compute Vh and Vl
            aH(end) = 1;
            aL(end) = 0;
            vH{end} = ((1-aH(end))*bisectionData{9,end} + aH(end)*bisectionData{8,end});
            vL{end} = ((1-aL(end))*bisectionData{9,end} + aL(end)*bisectionData{8,end});
            timeFull{end} = tEnd;
            
            for i = length(bisectionData)-1:-1:2
                
                linTime = bisectionData{12,i};
                alphaInd = bisectionData{7,i+1};
                Vhi = bisectionData{8,i};
                Vli = bisectionData{9,i};
                Vhi = Vhi(:,1:length(linTime))';
                Vli = Vli(:,1:length(linTime))';
                
                aH(i)     = aH(i+1)*alphaRef(alphaInd(2)) + (1-aH(i+1))*alphaRef(alphaInd(1));
                aL(i)     = aL(i+1)*alphaRef(alphaInd(2)) + (1-aL(i+1))*alphaRef(alphaInd(1));
                
                vH{i}    = ((1-aH(i))*Vli + aH(i)*Vhi)';
                vL{i}    = ((1-aL(i))*Vli  + aL(i)*Vhi)';
                timeFull{i} = linTime;
            end
            
            %only save the data which is less than teola - no extrapolation
            
            ind = find(tEnd2 <= teola);
            
            time{end} = tEnd2(ind);
            Vmeta{end} = (1-alphaM(end))*VLlast(:,ind) + alphaM(end)*VHlast(:,ind);
            
            
            for i = length(bisectionData)-3:-1:2
                
                linTime = bisectionData{12,i};
                alphaInd = bisectionData{7,i+1};               
                Vhi = bisectionData{8,i};
                Vli = bisectionData{9,i};
                Vhi = Vhi(:,1:length(linTime))';
                Vli = Vli(:,1:length(linTime))';
                alphaM(i) = alphaM(i+1)*alphaRef(alphaInd(2)) + (1-alphaM(i+1))*alphaRef(alphaInd(1));
                Vmeta{i} = ((1-alphaM(i))*Vli + alphaM(i)*Vhi)';      
                time{i} = linTime;
                
            end
            
            
            t = [time{:}];
            Vm = [Vmeta{:}];
            
            [time,ind,~] = unique(t);
            
            VmSpl = pchip(time,Vm(:,ind));
            
            VmSpl_t = @(t) ppval(VmSpl,t);
            
            t = [timeFull{:}];
            
            [timeFull,ind,~] = unique(t);
            Vh = [vH{:}];
            Vl = [vL{:}];
            VhSpl = pchip(timeFull,Vh(:,ind));
            VlSpl = pchip(timeFull,Vl(:,ind));
            
            VhSpl_t = @(t) ppval(VhSpl,t);
            VlSpl_t = @(t) ppval(VlSpl,t);
            
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Functions for coho
        
        % This function computes LDI models for input signals specified by Brockett annuli
        % region: current reachable region (an object of shape class)
        % states: the current region (1-4) for Brockett annulus of each signal.
        % vars:   index of interested nodes
        function [c,err] = dV_ldi(this,region,states,vars,varargin)
            if(nargin<4||isempty(vars)), vars = 1:this.synchronizerCCT.nodeNum; end
            vars = reshape(vars,1,[]);
            [c,err] = this.synchronizerCCT.dV_ldi(region,vars,varargin{:});
            % Override rather sum for input signals
            for i=1:length(this.sources)
                ind = this.map(i);
                if(any(vars==ind))
                    ibnds = region.project(ind).rec;
                    [cc,ee] = this.sources{i}.dV_ldi(ibnds,states(i),varargin{:});
                    mm = find(vars==ind); % index of vars
                    c(mm,[ind;end]) = repmat(cc,length(mm),1); err(mm) = repmat(ee,length(mm),1);
                end
            end
        end
        
        % The number of dimension of ODE
        function n = dim(this)
            n = this.synchronizerCCT.nodeNum;
        end
        % The name of the circuit
        function str = name(this)
            str = class(this.synchronizerCCT);
        end
        % design under test (circuit)
        function c = dut(this)
            c = this.synchronizerCCT;
        end
        % The number of input signals
        function n = inputNum(this)
            n = length(this.sources);
        end
        % voltage source of input signals
        function i = inputs(this,ind)
            if(nargin<2||isempty(ind))
                i = this.sources;
            else
                i = this.sources{ind};
            end
        end
        
        % Index of input signals
        function ind = inputNodes(this)
            ind = this.map;
        end
        % Index of circuit nodes
        % A wrapper of the circuit.find_port_index functions
        function inds = find_port_index(this,ports)
            inds = this.circuit.find_port_index(ports);
        end
        
        function model = getModel(this)
            model = this.defaultModel;
        end
        
        function vfull = getVfull(this,tSample,dataTransition,state)
            vfull = this.vfull(tSample,dataTransition,state);
        end
        
        function [Jac,dhda,vdot,dIdx_allDev,C_allDev] = computeJacVdot(this,t,dataTransition,metaStableState)
            %return the state with the data input held constant: vfull and vfullh
            %where the data input is driven from the clock source. This is to
            %seperate the homogenous and non-homogenous solutions to the circuit.
            
            [vdot,C,Itot,CapParams] = this.computeVdot(t,dataTransition,metaStableState);
            
            %AD is required for this analysis.
            if(~this.tbOptions.isAD)
                error('To continue this analysis, the intLAB package from Dr. Siegfried M. Rump found at: http://www.ti3.tu-harburg.de/rump/intlab/')
            end
            
            
            %now that the currents for each device is computed, we now have all
            %the information to compute the Jacobian, or partial Ids partial
            %Voltage nodes
            
            [Jac,dhda,dIdx_allDev,C_allDev] = this.computeJac(t,dataTransition,metaStableState,C,Itot,CapParams);
            
            dIdx_allDev = [dIdx_allDev{1};dIdx_allDev{2}];
            C_allDev    = [C_allDev{1};C_allDev{2}];
            
        end
        
        function [vdot,C,Idev,CapParams] = computeVdot(this,t,dataTransition,metaStableState)
            vfull = this.vfull(t,dataTransition,metaStableState);
            synchronizerCCTMapping = this.synchronizerCCTMap;
            
            
            %%%%%%
            %The circuit map returns a cell array where each entry contains is of
            %the form {ElementProperties,Map}, the elements are a cell array of all the
            %same kind of elements and the Map is the map that corresponds to the
            %appropriate mapping for those elements.
            %%%%%
            [~,numUniqueElements] = size(synchronizerCCTMapping);
            
            Itot = cell(1,numUniqueElements);
            Idev = cell(1,numUniqueElements);
            Ctot = cell(1,numUniqueElements);
            CapParams = cell(1,numUniqueElements);
            Ms = this.mapMatrix;
            for i = 1:numUniqueElements
                %get the circuit element grouping
                elementGroup = synchronizerCCTMapping{i};
                %extract the device properties and associated map
                deviceProperties = elementGroup{1};
                deviceMap = elementGroup{2};
                %extract the device type and model along with any specific model
                %parameters.
                deviceType = deviceProperties.deviceType;
                deviceModel = deviceProperties.deviceModel;
                %   deviceModelParams = deviceProperties.deviceModelParams;
                %struct which contains relevant parameters for the device in a
                %one to one correspondance. I.e.:
                %the parameter with field w {width_device1, width_device2, ...}
                vectorizedParams = deviceProperties.parameters;
                vectorizedParams.deviceMap = deviceMap;
                %determine the number of terminals per device and the number of
                %devices. Figure out the number of points to evaluate and shape V
                %such that each column represents a device for a particular
                %point.
                [numDeviceTerminals,numDevices] = size(deviceMap);
                [numCircuitNodes,numPointsToEval] = size(vfull);
                V = vfull(deviceMap,:);
                
                % V(which_terminal,which_point_in_phase_space,which_device)
                V = reshape(V,[numDeviceTerminals,numPointsToEval,numDevices]);
                %compute dvDevice in the form:
                %(numCircuitNodes,numPointsToEval,numDevices)
                %based on the user parameters, the input voltage: V and
                %simulation options. If the default model is not selected, the
                %model is extracted from the modelDictionary and then used to
                %compute dvDevice. This allows for multiple models to be used at
                %once during a simulation to support portions of circuits which
                %may have different model behaviors or other.
                if(strcmp(deviceModel,'default'))
                    [IDevice,capParams] = this.defaultModel.I(deviceType,vectorizedParams,V,this.tbOptions);
                    vectorizedParams.capParams = capParams;
                    CDevice = this.defaultModel.C(deviceType,vectorizedParams,V,this.tbOptions);
                else
                    model = this.modelDictionary.getModel(deviceModel,deviceModelParams);
                    [IDevice,capParams] = model.I(deviceType,vectorizedParams,V,this.tbOptions);
                    vectorizedParams.capParams = capParams;
                    CDevice = model.C(deviceType,vectorizedParams,V,this.tbOptions);
                end
                CapParams{i} = capParams;
                Idev{i} = IDevice;
                Itot{i} = Ms{i}*reshape(IDevice,numDeviceTerminals*numDevices*numPointsToEval,[]);
                Ctot{i} = Ms{i}*reshape(CDevice,numDeviceTerminals*numDevices*numPointsToEval,[]);
            end
            
            I = zeros(size(Itot{1}));
            C = zeros(size(Ctot{1}));
            
            for i = 1:length(Itot)
                C = C + Ctot{i};
                I = I + Itot{i};
            end
            
            %compute currents for voltage sources if any exist and their model
            %follows the default
            for i = 1:numUniqueElements
                %get the circuit element grouping
                elementGroup = synchronizerCCTMapping{i};
                %extract the device properties and associated map
                deviceProperties = elementGroup{1};
                deviceMap = elementGroup{2};
                %extract the device type and model along with any specific model
                %parameters.
                deviceType = deviceProperties.deviceType;
                deviceModel = deviceProperties.deviceModel;
                %make the currents and capacitance on each side of the voltage
                %source the same.
                if(strcmp(deviceType,'vsrc') && strcmp(deviceModel,'default'))
                    C(deviceMap,:) = C(deviceMap,:) + C(flipud(deviceMap),:);
                    I(deviceMap,:) = I(deviceMap,:) + I(flipud(deviceMap),:);
                end
            end
            
            isSource = this.is_src;

            I = I(~isSource);
            
            if(strcmp(this.tbOptions.capModel,'gnd'))
                C = C(~isSource);
                vdot = I./C;
            else
                
                %C here is the diagonal matrices of the sparse matrix I want to have,
                %MATLAB doesn't offer a convienient way to create a block
                %diagonal sparse matrix other than by the following:
                B = mat2cell(sparse(C),numCircuitNodes,numCircuitNodes);
                C = blkdiag(B{:});
                C = C(~isSource,~isSource);
                vdot = C\I;
            end
           
        end
        
        function [Jac,dhda,dIdxDev_Tot,Cdev_Tot] = computeJac(this,t,dataTransition,metaStableState,C,I,CapParams)
            
            
            vfull = this.vfull(t,dataTransition,metaStableState);
            synchronizerCCTMapping = this.synchronizerCCTMap;
            [~,numUniqueElements] = size(synchronizerCCTMapping);
            IResTot_f = cell(1,numUniqueElements);
            dIdx_tot = cell(1,numUniqueElements);
            Ih_tot = cell(1,numUniqueElements);
            dIdxDev_Tot = cell(1,numUniqueElements);
            Cdev_Tot = cell(1,numUniqueElements);
            Ms   = this.mapMatrix;
            %to do so we change the default model to use the AD version of the
            %model
            if(strcmp(this.tbOptions.capModel,'full'))
                Ctot = cell(1,numUniqueElements);
            end
            
            %%% Setup the voltage vector                
            vfullh = vfull;
            vfullh(this.inputSourceIndex,:) = dataTransition;
            

            Vf_AD = gradientinit(vfull);
            Vh_AD = gradientinit(vfullh);
            
            %set the voltage of the data transition node based on the time
            %and previously determined data transition time
            
            din = this.inputSource.V(t,Vh_AD(this.inputSourceIndex,:));
            Vh_AD(this.inputSourceIndex,:) = din;
            
            for i = 1:numUniqueElements
                %get the circuit element grouping
                elementGroup = synchronizerCCTMapping{i};
                %extract the device properties and associated map
                deviceProperties = elementGroup{1};
                deviceMap = elementGroup{2};
                %extract the device type and model along with any specific model
                %parameters.
                deviceType = deviceProperties.deviceType;
                deviceModel = deviceProperties.deviceModel;
                %   deviceModelParams = deviceProperties.deviceModelParams;
                %struct which contains relevant parameters for the device in a
                %one to one correspondance. I.e.:
                %the parameter with field w {width_device1, width_device2, ...}
                vectorizedParams = deviceProperties.parameters;
                vectorizedParams.deviceMap = deviceMap;
                %determine the number of terminals per device and the number of
                %devices. Figure out the number of points to evaluate and shape V
                %such that each column represents a device for a particular
                %point.
                [numDeviceTerminals,numDevices] = size(deviceMap);
                [numCircuitNodes,numPointsToEval] = size(vfull);
                
                %extract what was computed from the transistor current
                %evaluator, the current multiplication factor is [-1;0;1;0]
                %where:
                % factor 1: Current entering/leaving Drain terminal
                % factor 2: Current entering/leaving Gate terminal
                % factor 3: Current entering/leaving the Source terminal
                % factor 4: Current entering/leaving the Body terminal
                Idev_f = I{i};
                IdsDev_f = reshape(Idev_f,[numDeviceTerminals,numPointsToEval,numDevices]);
                sign_f = sign(IdsDev_f);
                sign_f = sign_f(3,:,:);
                sign_f = reshape(sign_f,1,[]);
                IdsDev_f = abs(IdsDev_f(3,:,:));
                IdsDev_f = reshape(IdsDev_f,1,[]);
                
                Idev_h = I{i};
                IdsDev_h = reshape(Idev_h,[numDeviceTerminals,numPointsToEval,numDevices]);
                IdsDev_h = abs(IdsDev_h(3,:,:));
                IdsDev_h = reshape(IdsDev_h,1,[]);

                
                V = Vf_AD(deviceMap,:);
                Vh = Vh_AD(deviceMap,:);
                
                % V(which_terminal,which_point_in_phase_space,which_device)
                V = reshape(V,[numDeviceTerminals,numPointsToEval,numDevices]);
                Vh = reshape(Vh,[numDeviceTerminals,numPointsToEval,numDevices]);
                
                if(strcmp(deviceModel,'default'))
                    [IdevRes_f,dIdx] = this.defaultModelAD.IdevRes(deviceType,vectorizedParams,V,[sign_f;IdsDev_f],numCircuitNodes,this.tbOptions);
                    Idev_h = this.defaultModelAD.I(deviceType,vectorizedParams,{Vh,IdsDev_h'},this.tbOptions);
                    if(strcmp(this.tbOptions.capModel,'full'))
                        vectorizedParams.capParams = CapParams{i};
                        Cnew = this.defaultModelAD.C(deviceType,vectorizedParams,Vh,this.tbOptions);
                    end
                    
                else
                    model = this.modelDictionary.getModel(deviceModel,deviceModelParams);
                    [IdevRes_f,dIdx] = model.IdevRes(deviceType,vectorizedParams,V,[sign_f;IdsDev_fAD],numCircuitNodes,this.tbOptions);
                    Idev_h = model.I(deviceType,vectorizedParams,{Vh,IdsDev_hAD},this.tbOptions);
                    if(strcmp(this.tbOptions.capModel,'full'))
                        vectorizedParams.capParams = CapParams{i};
                        Cnew = model.C(deviceType,vectorizedParams,Vh,this.tbOptions);
                    end
                end
                
                
                IResTot_f{i} = Ms{i}*reshape(IdevRes_f,numDeviceTerminals*numDevices*numPointsToEval,[]);
               
                Ih_tot{i} = Ms{i}*reshape(Idev_h,numDeviceTerminals*numDevices*numPointsToEval,[]);
                
                dIdx_tot{i} = Ms{i}*dIdx;
                dIdxDev_Tot{i} = dIdx;
                Cdev_Tot{i} = C;
                if(strcmp(this.tbOptions.capModel,'full'))
                    Ctot{i} = Ms{i}*Cnew;
                    Cdev_Tot{i} = Cnew.x;
                end

            end
            I_h = zeros(size(Ih_tot{1}));
            Ires_f = zeros(size(IResTot_f{1}));
            dIdx = zeros(size(dIdx_tot{1}));
            if(strcmp(this.tbOptions.capModel,'full'))
                C = zeros(size(Ctot{1}));
            end
            
            for i = 1:length(Ih_tot)
                I_h = I_h + Ih_tot{i};
                Ires_f = Ires_f + IResTot_f{i};
                dIdx    = dIdx + dIdx_tot{i};
                if(strcmp(this.tbOptions.capModel,'full'))
                    C = C + Ctot{i};
                end
            end
            
            %compute currents for voltage sources if any exist and their model
            %follows the default
            for i = 1:numUniqueElements
                %get the circuit element grouping
                elementGroup = synchronizerCCTMapping{i};
                %extract the device properties and associated map
                deviceProperties = elementGroup{1};
                deviceMap = elementGroup{2};
                %extract the device type and model along with any specific model
                %parameters.
                deviceType = deviceProperties.deviceType;
                deviceModel = deviceProperties.deviceModel;
                %make the currents and capacitance on each side of the voltage
                %source the same.
                if(strcmp(deviceType,'vsrc') && strcmp(deviceModel,'default'))
                    Ires_f(deviceMap,:) = Ires_f(deviceMap,:) + Ires_f(flipud(deviceMap),:);
                    I_h(deviceMap,:) = I_h(deviceMap,:) + I_h(flipud(deviceMap),:);
                    C(deviceMap,:)   = C(deviceMap,:) + C(flipud(deviceMap),:);
                end
            end
            isSource = this.is_src;
            I_h = I_h(~isSource);
            dIdx   = dIdx(~isSource,~isSource);
            if(strcmp(this.tbOptions.capModel,'gnd'))
                h = I_h./C;
                Jac = dIdx./C';
                dhda = h.dx(:,this.inputSourceIndex);
            else
                Cap = C.x;
                Cap = Cap(~isSource,~isSource);
                f = Cap\I_h.x;
                dC = C.dx;
                Idx = I_h.dx;
                
                %compute partial f / partial alpha
                dCda = dC(~isSource,~isSource,this.inputSourceIndex);
                dhda = Cap\(dCda*f+Idx(:,this.inputSourceIndex));
                
                dCdx = zeros(length(f),length(f));
                offset = length(this.sources);
                    for i = 1:length(f)
                        dCdx(i,:) = dC(~isSource,~isSource,i+offset)*f;
                    end
                
                Jac = Cap\(dCdx + dIdx);
            end
        end
        
    function maps = computeMapMatrix(this)
        
        devMaps = this.synchronizerCCTMap;
        maps = cell(1,length(devMaps));
        maxNodeNum = 0;
        for i = 1:length(devMaps)
            m = reshape(devMaps{i}{2},1,[]);
            if(max(m) > maxNodeNum)
                maxNodeNum = max(m);
            end
        end
        
        for i = 1:length(devMaps)
            m = reshape(devMaps{i}{2},1,[]);
            m = m';
            Mdev = sparse(m,1:size(m),ones(size(m)));
            [x,~] = size(Mdev);
            % this is a hacky solution which is not fast, BUT this mapping
            % only gets computed once
            if(x < maxNodeNum)
                Mdev(maxNodeNum,:) = 0;
            end
            maps{i} = Mdev;
        end
    end
        
        function Mdev = computeMapMatrixDevice(this,in,out)
            Mtotal = [this.mapMatrix{:}];
            devMaps = this.synchronizerCCTMap;
            maps = cell(1,length(devMaps));
            for i = 1:length(devMaps)
                devSet_i = devMaps{i}{2};
                devMapTemp = zeros(size(devSet_i));
                for j = 1:length(devSet_i)
                    if(any(ismember(devSet_i(:,j),in)) && any(ismember(devSet_i(:,j),out)))
                        devMapTemp(:,j) = 1;
                    end
                end
                maps{i} = reshape(devMapTemp,1,[]);
                
            end
            
            mapr = [maps{:}];
            mapr = logical(mapr);
            Mdev = zeros(size(Mtotal));
            Mdev(:,mapr) = Mtotal(:,mapr);
            Mdev = sparse(Mdev);
        end
        
        function Mdev = computeMapMatrixDeviceTx(this,in,out,type)
            Mtotal = [this.mapMatrix{:}];
            devMaps = this.synchronizerCCTMap;
            maps = cell(1,length(devMaps));
            for i = 1:length(devMaps)
                deviceType = devMaps{i}{1}.deviceType;
                devSet_i = devMaps{i}{2};
                devMapTemp = zeros(size(devSet_i));
                if(strcmp(deviceType,type))
                    for j = 1:length(devSet_i)
                        if(any(ismember(devSet_i(:,j),in)) && any(ismember(devSet_i(:,j),out)))
                            devMapTemp(:,j) = 1;
                        end
                    end
                end
                maps{i} = reshape(devMapTemp,1,[]);
                
            end
            
            mapr = [maps{:}];
            mapr = logical(mapr);
            Mdev = zeros(size(Mtotal));
            Mdev(:,mapr) = Mtotal(:,mapr);
            Mdev = sparse(Mdev);
        end
        
    end
    
    methods(Access='private')
        % compute the voltage for all circuit nodes
        % v_ode: each column for a time
        function vfull = vfull(this,t,dataTransition,v_ode)
            assert(length(t)==size(v_ode,2));
            
            v_odeShape = size(v_ode);
            
            % v_ode = reshape(v_ode,[],this.tbOptions.numParallelCCTs);
            % voltages of sources at time t
            vs = this.eval_vs(t);%./this.tbOptions.capScale);
            % assign voltage source
            iss = this.is_src;
            vfull = zeros(this.synchronizerCCT.nodeNum,v_odeShape(2)); %the full vector for one circuit
            vfull(this.map,:) = vs; % fill in voltage sources map to circuit port
            %vfull = repmat(vfull,this.tbOptions.numParallelCCTs,1);   %get num parallel cct copies
            %vfull = reshape(vfull,[],this.tbOptions.numParallelCCTs);
            vfull(~iss,:) = v_ode;
            
            synchronizerMapping = this.synchronizerCCTMap;
            %%%%%%
            %The circuit map returns a cell array where each entry contains is of
            %the form {ElementProperties,Map}, the elements are a cell array of all the
            %same kind of elements and the Map is the map that corresponds to the
            %appropriate mapping for those elements.
            %%%%%
            [~,numUniqueElements] = size(synchronizerMapping);
            
            for i = 1:numUniqueElements
                dev = synchronizerMapping{i};
                devParams = dev{1};
                devMap = dev{2};
                if(strcmp(devParams.deviceType, 'vsrc'))
                    internalVsrc = dev{3};
                    for j = 1:length(internalVsrc)
                        src = internalVsrc{j};
                        vSrcMap = devMap(:,j);
                        v = vfull(vSrcMap,:);
                        vfull(vSrcMap(1),:) = v(2,:) + src.V(t);%./this.tbOptions.capScale);
                    end
                end
            end
            
            vfull(this.inputSourceIndex,:) = this.inputSource.V(t,dataTransition);
            
        end
        % compute voltage for circuit inputs
        % vs: each column for a time
        function vs = eval_vs(this,t)
            vs = zeros(length(this.sources),length(t));
            for i=1:length(this.sources)
                vs(i,:) = this.sources{i}.V(t);
            end
        end
        % index of source nodes
        function iss = is_src(this)
            n = this.synchronizerCCT.nodeNum;
            iss = false(n,1);
            iss(this.map) = true;
        end
        
        function setInterval(this,newInterval)
            this.dinInterval = newInterval;
        end
             
        function dwdt = dwdt_ode(~,jac,w_ode)
            dwdt = -w_ode'*jac;
            dwdt = dwdt';
        end
           
        
        function clockEdgeTimes = findClockEdges(this,t)
            
            clkEdgeOptions = this.tbOptions.clockEdgeSettings;
            
            %produce a clock trace and find where the clock edge crosses
            %between 'low' percentage and 'high' percentage of vdd.
            clk = this.clockSource.V(t);
            edgeIndex = clk > this.tbOptions.vdd*clkEdgeOptions.low...
                & clk < this.tbOptions.vdd*clkEdgeOptions.high;
            
            edges = t(edgeIndex);
            [~,iter] = size(edges);
            %this is to ensure that the last edge is captured.
            if(mod(iter,2)==0)
                iter = iter+1;
                edges(iter)= 0;
            end
            
            edgeNum = 1;
            counter = 1;
            timeSum = 0;
            
            clockEdgeTimes = zeros(1,clkEdgeOptions.numberOfEdges);
            
            for i = 1:(iter-1)
                %if the differences between the edges is less than clkEdgeSpread apart
                %then increase the timeSum and counter otherwise they belong to
                %a new edge
                if(abs(edges(i)-edges(i+1)) < clkEdgeOptions.clkEdgeSpread)
                    timeSum = edges(i) + timeSum;
                    counter = counter + 1;
                else
                    clockEdgeTimes(edgeNum) = (timeSum + edges(i))/counter;
                    %if the number of edges requested is met, exit. Otherwise
                    %reset the counter and timeSum
                    if(edgeNum > clkEdgeOptions.numberOfEdges)
                        break;
                    end
                    edgeNum = edgeNum + 1;
                    counter = 1;
                    timeSum = 0;
                end
            end
        end
        
        function [g,dgdt,lambda] = gainSync(~,wNorm,jac,beta,dhda)
            lambda = wNorm'*jac*wNorm;
            g = wNorm'*beta;
            dgdt = wNorm'*dhda+lambda*g;
            
        end
        
        function [value,isterminal,direction] = betaRenorm(this,t,beta)
            betaOptions = this.tbOptions.betaSimOptions;
            
            isterminal = 1;
            direction = 1;
            
            if(t < betaOptions.minSimTime)
                a = -1;
            end
            
            value = norm(beta) - betaOptions.threshold;
        end
        
    end
end
