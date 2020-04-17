% testbench is a superclass of nestedBisection, nestedBisection has the
% following properties:
%   inputSource: the input to the synchronizer, must implement V(t,tin),
%   where tin is a parameter which determines when the input makes a
%   transistion from low to high (or high to low).
%   dinInterval: The interval within which the data transistion is supposed
%   to occur in to continue the nested bisection algorithm
%   inputSourceIndex: The index in the circuit description which
%   corresponds to the input source
%   clockSource: The clock source for the synchronizer
%   nodeMask: The set of nodes which are used to bisect on corresponding to
%   the rising/falling edges of the clock.
%
% constructor
%   nestedBisection(synchronizerCCT,ports,sources,dinInterval,dinSource,...
%                   clockSource,nodeMask,uVcrit,modelDictionary,...
%                   defaultModel,defaultModelProcess,varargin)
%       synchronizerCCT: The synchronizer circuit definition
%       ports:          a cell of "node" objects corresponding to source
%                       nodes.
%       sources:        a cell of "circuit" objects which must implement
%                       their V(t) function.
%       dinInterval:    The initial interval for tin where one end point
%                       must correspond to the synchronizer to output a low
%                       signal and the other end point correspond to the
%                       synchronizer to output a high signal
%       dinSource:      The input source whihc must implement a
%                       parameterized source voltage V(t,tin)
%       clockSource:    The source used to clock the synchronizer circuit
%       nodeMask:       A list of nodes which correspond to which node is
%                       used to determine when the synchronzier has settled
%                       high and low.
%   inputSource: the input source which must implement V(t,tin) where tin
%     is the time of the input transistion to the synchronizer circuit
%     din
%     ports:    an cell of "node" objects
%     sources: an cell of "circuit" objects which must implement the V(t) function.
%          the size of ports and sources must be the same.
%          the voltage source is applied to the corresponding port
% methods:
%  bisectionRuns = simulate(filename,tspan,tCrit,v0)
%       This function launches the nested bisection algorithm with the
%       apporpriate settings defined in tbOptions.
%
%           filename: The filename given where the resultant trajectories
%           and diagnostic data is saved.
%
%           tspan: The time span of the simulation time, t=0 is shifted
%           to be at the time at which the first clock edge occurs as
%           defined by the clock source of the synchronizer.
%
%           tCrit: The user specified time by which the synchronizer is
%           supposed to be at well defined logical values. When the
%           simulation time exceeds tCrit, the nested bisection algorithm
%           checks if the user defined stopping conditions have been met.
%
%           OPTIONAL SETTINGS:
%               v0: A user specified initial condition to apply to all
%               simultaneous synhcronizers. Default is to set all nodes to
%               0.
%
%           bisectionRuns: The data structure returned by the nested
%           bisection algorithm which contains the following results in a
%           cell array:
%               bisectionRuns{1,:}: [tinH,tinL] for each epoch
%               bisectionRuns{2,:}: the saved fraction by which the
%               algorithm shrank the interval at each epoch
%               bisectionRuns{3,:}: [startTime,eventTime] the start and end
%               time of each epoch
%               bisectionRuns{4,:}: [Vh,Vl] the circuit state for the
%               trajectory which settled high and low at the end of each
%               epoch
%               bisectionRuns{5,:}: The epoch number
%               bisectionRuns{6,:}: A centered numerical difference
%                                   estimate for beta(t)
%               bisectionRuns{7,:}: [indL,indH] the trajectory index which
%               corresponds to the trajectory which settled low and the one
%               which settle high for the given epoch.
%               bisectionRuns{8,:}: The trajectory for the circuit which
%               settled high between the start time and the end for the
%               given epoch.
%               bisectionRuns{9,:}: The trajectory for the circuit which
%               settled low between the start time and the end for the
%               given epoch.
%               bisectionRuns{10,:}: The time points for the given epoch.
%               bisectionRuns{11,:}: The difference between tinH and tinL
%               for the given epoch.
%               bisectionRuns{12,:}: The time interval where the circuit
%               has linear behaviour during the given epoch.
%               bisectionRuns{13,:}: Saved ODE points for each epoch if
%               debug mode is enabled.
%               bisectionRuns{14,:}: The saved voltage states for each
%               circuit on the given epoch.
%
%   dv = dV_ode(t,v_ode): construct ODE for ODE nodes
%     t:      used to evaluates voltage of sources
%     v_ode:  voltage of ODE nodes
%     dv_ode: time derivative of ODE voltages
%
%   clockEgeTimes = findClockEdges(t): returns when the clock signal is
%   within the specified window defined in tbOptions
%
%     t: the input range to search for clock edges (rising or falling).
%
%   private methods:
%
%   index = findNextBracket(Vin,t,edgeTimes): Function which returns the
%   indices of the circuit states upon which the next epoch is subdivided
%
%     Vin: The entire voltage state for all simulataneously simulated
%     circuits
%
%     t: The time for the given epoch to determine which circuit node is
%     considered for determining the indices
%
%     edgeTimes: Set of times at which the clock signal has reached the
%     desired time interval specified by tbOptions.
%
%   [startTime,newStates,VL,VH,diffStates,timeInt] =
%   findNextStartStatesAndTime(Vin,t,index): returns when the time at which
%   the linearity assumption set in tbOptions is violated, the new initial
%   conditions for the next epoch, the circuit state at the end of the
%   linearity assumption for both the states which eventually settle high
%   and low, and finally the time interval where the linearity assumption
%   holds.
%
%       Vin: The input voltage state of all simultaneously simulated
%       synchronizers
%
%       t: the time points for the particular simulation run
%
%       index: The computed indices to determine which circuit settles high
%       and which one settles low.
%
%       startTime: The new start time of the nested bisection epoch.
%
%       newStates: The new starting conditions for the next epoch.
%
%       VL: The end point at the end of the epoch which eventually settles
%       low
%
%       VH: The end point at the end of the epoch which eventually settles
%       high
%
%       diffStates: The difference between the true circuit states and the
%       assumed linear model.
%
%       timeInt: The interval of time where the linearity assumption holds.
%
%   [value,isterminal,direction] = stopCriteria(t,V): Stopping function for
%   the integrator. Makes a call to an external user provided function
%   stopToDigitalSim(t,Vin,vDot,digitalSimOptions) which determines if the
%   user specified criteria is met to stop the numerical integration.
%
%       t: The time computed by the integrator
%
%       V: The integrator voltage state
%
%       value: the value returned by the stopCriteria function
%       isterminal: determines if the condition is terminal
%       direction:  determines what direction the value needs to have
%       changed by
%
%       The above implements the event function from MATLABs ode
%       integrators
%
%
%   bisectionPlotV(plotHandles,t,Vin): This function plots the trajectories
%   of the nodes provided by the user. Each subplot displays the
%   trajectories for all simultaneously simulated synchronizer circuits.
%
%       plotHandles: the plot handles created at setup time
%
%       t: The time points simulated from numerical integration
%
%       Vin: The voltage states of the synchronizer
%
%   bisectionPlotBeta(plotHandles,tbeta,beta): This function plots the
%   numercial approximation of beta(t) at each nested bisection epoch
%
%       plotHandles: the plot handles createdat setup time
%
%       tbeta: the time for which the synchronizer is behaving in linear
%       fashion, tbeta is determined by findNextStartStatesAndTime (see
%       above)
%
%       beta: the associated numerical estimate for beta
%
%   plotHandles = bisectionPlotSetup(tspan,tCrit,tOff): reads in any
%   specific plotting options defined by tbOptions.plotOptions. By default
%   there are no options and this function creates a plot where each node
%   defined in the mask is plotted on its own subplot and creates a plot to
%   show the numerical estimation of beta. A dotted line shows the user's
%   critical time selection and tOff time shifts the plot so that t=0
%   occurs at the first clock edge.
%
%       tspan: the simulation time span
%
%       tCrit: the user selected critical time by which time the
%       synchronizer needs to have settled.
%
%       tOff: the computed offset to line up t=0 to be the time at which
%       the first clock edge occurs.
%
%       OPTIONAL: Provisions are supplied but not currently implemented to
%       add additional plotting options.
%
%   vfull = vfull(t,v_ode): A function which adds extenal voltage sources
%   to the computed integrated circuit state. The full voltage vector is
%   necessary for properly biasing all transistor devices in the circuit.
%
%       t: The current integrated time step
%
%       v_ode: The current integrated voltage state of the circuit state
%
%   setInterval(newInterval): A function which sets the time interval of
%   the input transistion time for the synchronizer's input source.
%
%       newInterval: The interval for [tHi,tLow] for the next bisection
%       run.

classdef nestedBisection < testbench
    properties (GetAccess='public', SetAccess='private')
        
        %the input source which will be used as a variable for bisection even
        %though bisection is done on voltage state the input source needs a
        %variable with which the delay can vary, a requirement to the method.
        inputSource = [];
        dinInterval = [];
        inputSourceIndex = 0;
        clockSource = [];
        
        %nodeMask is an array of node indices that correspond to nodes to check
        %in determining the outcome of the synchronizer. i.e. nodeMask(1) is
        %the node to look at during the first clock edge to determe if the
        %synchronizer will settle low or high after N clock edges
        nodeMask = [];
        
        
    end
    properties (GetAccess = 'protected', SetAccess = 'private')
        synchronizerCCTMap;
        %total circuit mapping
        syncMatrixMap = [];
    end
    
    methods
        % to create a nestedBisection testbench you need to provide:
        % - circuit: a synchronizer circuit description
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
        %                                  - 'vdd': the circuit's overall
        %                                  vdd voltage level.
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
        
        function this = nestedBisection(synchronizerCCT,ports,sources,dinInterval,dinSource,...
                clockSource,nodeMask,modelDictionary,defaultModel,defaultModelProcess,varargin)
            if(nargin<1), error('not enough parameters'); end
            if(nargin<2), ports = {}; end
            if(nargin<3), sources = {}; end
            if(nargin<4), error('a model dictionary is required'); end
            if(nargin<10), defaultModel = 'MVS'; defaultModelProcess = 'PTM 45nmHP'; end
            if(nargin<11)
                integratorOptions = odeset('RelTol',1e-6,'AbsTol',1e-6);
                
                %%%%%% the following default options are for nested bisection of
                %%%%%% the synchronizer analysed in the ASYNC 2018 paper, i.e. a
                %%%%%% 2-Flop passgate synchronizer. These options for your
                %%%%%% synchronizer should be provided.
                digitalSimOptions = struct('vdd',1.0,'minSimTime',6e-10,'threshold',0.1,...
                    'stopGain',30);
                
                clkEdgeSettings = struct('low',0.48,'high',0.52,'clkEdgeSpread',20e-12,'numberOfEdges',3);
                
                tbOptions = struct('capModel','full','capScale',1e10,'vdd',1.0,'temp',298,'numParallelCCTs',10,...
                    'stepOut',1,'polyFitDegree',3,'polyFitError',1e-1,'isAD',false,'transSettle',1e-10,...
                    'plotOptions',false,'integratorOptions',integratorOptions,'clockEdgeSettings',clkEdgeSettings,...
                    'digitalSimOptions',digitalSimOptions,'integrator',@ode45,...
                    'debug',false);
            end
            
            if(length(ports)~=length(sources))
                error('incorrect size of sources and ports');
            end
            
            if(isempty(tbOptions))
                tbOptions = varargin;
            end
            
            this = this@testbench(synchronizerCCT,ports,sources,modelDictionary,defaultModel,defaultModelProcess,tbOptions);
            tbOptions = this.getTestbenchOptions;
            tbOptions.integratorOptions.Events = @this.stopCriteria;
            
            %this just makes it so that the data transition time which happens
            %first causes the synchronizer to settle low.
            if(dinInterval(1) > dinInterval(2))
                this.dinInterval = dinInterval*this.tbOptions.capScale;
            else
                this.dinInterval = fliplr(dinInterval)*this.tbOptions.capScale;
            end
            this.inputSource = dinSource{2};
            this.clockSource = clockSource;
            %subtract the source nodes from the mask.
            this.nodeMask = nodeMask - length(sources);
            this.synchronizerCCTMap = this.circuit.getCircuitMap;
            if(any(this.map==0))
                error('can not find ports');
            end
            
            %find the index at which the input source is evaluated at
            this.inputSourceIndex = this.circuit.find_port_index(dinSource{1});
            %take the uVcrit vector and restrict it to state variables only.
            %      tbOptions.tCritSettings.uVcrit = tbOptions.tCritSettings.uVcrit - length(this.sources);
            
            this.setTestbenchOptions(tbOptions);
            
            % computes the matrix map to simulate the number of specified
            % parallel circuits
            maps = this.computeMapMatrix;
            M = cell(1,length(maps));
            for i = 1:length(maps)
                M{i} = sparse(kron(eye(this.tbOptions.numParallelCCTs),maps{i}));
            end
            this.syncMatrixMap = M;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Functions for simulation
        
        % This function simulate the circuit
        function [bisectionRuns] = simulate(this,filename,tspan,tCrit,v0)
            
            if(ischar(filename)||isstring(filename))
                fprintf('starting nested bisection and saving data to filename: %s ...\n',filename);
            else
                error('filename needs to be a string or array of characters')
            end
            n = this.circuit.nodeNum;
            tspan = tspan*this.tbOptions.capScale;
            tCrit = tCrit*this.tbOptions.capScale;
            % initial values
            if(nargin < 5)
                v0 =zeros(n,1); % initialize all nodes to zero
            elseif(iscell(v0)) % user specified a few nodes to initialize
                v00 = v0;
                v0 = zeros(n,1); % set the others to 0
                v0(v00{1}) = v00{2};
            end
            % Integrator
            %set the default integrator options if none are specified
            opts = this.tbOptions.integratorOptions;
            integrator = this.tbOptions.integrator;
            
            % simulate ODE nodes
            isOde = ~this.is_src;
            v0_ode = v0(isOde);
            v0_ode = repmat(v0_ode,[this.tbOptions.numParallelCCTs,1]); %set up n parallel circuits to run concurrently
            
            initialInterval = this.dinInterval;
            
            if(isempty(v0_ode))
                t = reshape(tspan, [], 1);
                v0_ode = zeros(length(t), 0);
            else
                %the following is to find where the clock edges are, split the
                %time interval into approximately 1ps resolution
                res = round((tspan(2)-tspan(1))/(1e-12*this.tbOptions.capScale));
                t = linspace(tspan(1),tspan(2),res);
                
                clockEdges = this.findClockEdges(t);
                startTime = tspan(1);
                stopTime = tspan(2);
                bisectionData = cell(14,100);
                %name the fields of saved data
                bisectionData{1,1} = '[alphaH,alphaL]';
                bisectionData{2,1} = 'fraction';
                bisectionData{3,1} = '[startTime,eventTime]';
                bisectionData{4,1} = '[Vh at eventTime, Vl at eventTime]';
                bisectionData{5,1} = 'bisection count';
                bisectionData{6,1} = 'beta numerical approx';
                bisectionData{7,1} = '[traceNumLow, traceNumHigh]';
                bisectionData{8,1} = 'trajectory high';
                bisectionData{9,1} = 'trajectory low';
                bisectionData{10,1} = 'time points';
                bisectionData{11,1} = 'delta alpha';
                bisectionData{12,1} = 'linear time interval';
                bisectionData{13,1} = 'ode points';
                bisectionData{14,1} = 'computed trajectories';
                count = 1;
                runningFrac = 1;
                clk = @(t) this.clockSource.V(t) - 0.5*this.tbOptions.vdd;
                tOff = fzero(clk,clockEdges(1)/this.tbOptions.capScale);
                numStatesPerCCT = length(v0_ode)/this.tbOptions.numParallelCCTs;
                plotHandles = this.bisectionPlotSetup(tspan./this.tbOptions.capScale,tCrit/this.tbOptions.capScale,tOff);
                
                tw = 1;
                
                stop = false;
                
                fprintf('Starting bisection...\n')
                tic()
                while(~stop)
                    
                    if(this.tbOptions.debug)
                        ind = 1;
                        data = cell(5,5000);
                        this.setOdeSave(ind,data);
                    end
                    
                    [t, v_ode,timeEvent,vEvent,indexEvent] = integrator(@(t, v_ode)(this.dV_ode(t, v_ode)), [startTime stopTime], v0_ode, opts);
                    
                    indAfterSettle = t > this.tbOptions.transSettle*this.tbOptions.capScale;
                    t = t(indAfterSettle);
                    tEnd = t(end);
                    v_ode = v_ode(indAfterSettle,:);
                    this.bisectionPlotV(plotHandles(1),t./this.tbOptions.capScale-tOff,v_ode)
                    
                    
                    index = this.findNextBracket(v_ode(end,:),t(end),clockEdges);
                    %this is to prevent state bisection to occur on transients
                    %due to initial conditions
                    
                    [startTime,v0_ode,VL,VH,diffStates,linTint] = this.findNextStartStatesAndTime(v_ode,t,index);
                    count = count +1;
                    alphaValues = linspace((initialInterval(1)/runningFrac),(initialInterval(2)/runningFrac),this.tbOptions.numParallelCCTs)./this.tbOptions.capScale;
                    deltaAlpha = alphaValues(1)-alphaValues(2);
                    runningFrac = abs(runningFrac*length(alphaValues)/(index(1)-index(2)));
                    tw = 1/runningFrac;
                    betaApprox = vecnorm(diffStates./(3*deltaAlpha));
                    
                    %save the trajectory which went high and the trajectory
                    %which went low
                    vTrajH = v_ode(:,1+numStatesPerCCT*(index(2)-1):numStatesPerCCT*index(2))';
                    vTrajL = v_ode(:,1+numStatesPerCCT*(index(1)-1):numStatesPerCCT*index(1))';
                    
                    %if the simulation ends at a time greater than the
                    %specified tCrit, call the function nestedBisectionStop to
                    %determine if the stopping condition has been met.
                    if(tEnd > tCrit)
                        stop = nestedBisectionStop(this.tbOptions,numStatesPerCCT,v_ode,tCrit,t);
                    end
                    
                    this.bisectionPlotBeta(plotHandles{2},linTint./this.tbOptions.capScale-tOff,betaApprox);
                    
                    t = t';
                    t = t./this.tbOptions.capScale;
                    linTint = linTint';
                    linTint = linTint./this.tbOptions.capScale;
                    bisectionData{1,count} = this.dinInterval./this.tbOptions.capScale;
                    bisectionData{2,count} = (index(2) - index(1))/(this.tbOptions.numParallelCCTs - 1);
                    bisectionData{3,count} = [startTime,timeEvent]./this.tbOptions.capScale;
                    bisectionData{4,count} = [VH,VL];
                    bisectionData{5,count} = count;
                    bisectionData{6,count} = diffStates./(deltaAlpha);
                    bisectionData{7,count} = index;
                    bisectionData{8,count} = vTrajH;
                    bisectionData{9,count} = vTrajL;
                    bisectionData{10,count} = t;
                    bisectionData{11,count} = deltaAlpha;
                    bisectionData{12,count} = linTint;
                    bisectionData{13,count} = this.getOdeSave;
                    bisectionData{14,count} = v_ode;
                    
                    
                end
                toc()
                fprintf('Bisection simulation completed with a relative tw = %d \n',tw)
                %save the data
                bisectionRuns = bisectionData(:,1:count);
                save(strcat(filename,'.mat'),'bisectionRuns');
                %reset the interval to return the test bench to the original
                %initial conditions
                this.setInterval(initialInterval./this.tbOptions.capScale);
                
            end
            
        end % simulate
        
        % This function is a wrapper of circuit.dV by setting input values
        % as source.V(t) for the Matlab integrator.
        function dv = dV_ode(this,t,v_ode)
            vodeShape = size(v_ode);
            vfull = this.vfull(t,v_ode);
            synchronizerCCTMapping = this.synchronizerCCTMap;
            
            %%%%%%
            %The circuit map returns a cell array where each entry contains is of
            %the form {ElementProperties,Map}, the elements are a cell array of all the
            %same kind of elements and the Map is the map that corresponds to the
            %appropriate mapping for those elements.
            %%%%%
            [~,numUniqueElements] = size(synchronizerCCTMapping);
            Itot = cell(1,numUniqueElements);
            Ctot = cell(1,numUniqueElements);
            Ms   = this.syncMatrixMap;
            
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
                    if(strcmp(deviceType,'vsrc'))
                        CDevice = CDevice(flipud(deviceMap),:);
                        IDevice = IDevice(flipud(deviceMap),:);
                    end
                else
                    model = this.modelDictionary.getModel(deviceModel,deviceModelParams);
                    IDevice = model.I(deviceType,vectorizedParams,V,this.tbOptions);
                    CDevice = model.C(deviceType,vectorizedParams,V,this.tbOptions);
                end
                
                
                %the following is unique to mos devices for the bisection
                %algorithm, the following flattens the currents back into the
                %form [I_cct1,I_cct2,...,I_cctN] where I_cct1 are all the
                %currents for the first copy of the circuit.
                
                % Example a 3 device circuit: I_cct1 = [Idev1,Idev2,Idev3] where
                % Idev1 = [Id;Ig;Is;Ib]
                
                %Store the result as [Idev_1,Idev_2,...Idev_n]
                %Similar applies for the Capacitance.
                Itot{i} = Ms{i}*reshape(IDevice,numDeviceTerminals*numDevices*numPointsToEval,[]);
                Ctot{i} = Ms{i}*reshape(CDevice,numDeviceTerminals*numDevices*numPointsToEval,[]);
            end
            
            I = zeros(size(Itot{1}));
            C = zeros(size(Ctot{1}));
            for i = 1:length(Itot)
                C = C + Ctot{i};
                I = I + Itot{i};
            end
            iss = this.is_src;
            isSource  = repmat(iss,this.tbOptions.numParallelCCTs,1);
            Inodes = I(~isSource);
            C =  this.tbOptions.capScale * C;
            
            if(strcmp(this.tbOptions.capModel,'gnd'))
                Cnodes = C(~isSource);
                dv = Inodes./Cnodes;
            else
                
                %C here is the diagonal matrices of the sparse matrix I want to have,
                %MATLAB doesn't offer a convienient way to create a block
                %diagonal sparse matrix other than by the following:
                B = mat2cell(sparse(C),numCircuitNodes*ones(1,this.tbOptions.numParallelCCTs),numCircuitNodes);
                C = blkdiag(B{:});
                Cnodes = C(~isSource,~isSource);
                dv = Cnodes\Inodes;
            end
            dv = reshape(dv,vodeShape);
            if(this.tbOptions.debug)
                odeSave = this.getOdeSave;
                ind = odeSave{1};
                data = odeSave{2};
                data{1,ind} = C(~iss,:);
                data{2,ind} = I(~iss,:);
                data{3,ind} = dv;
                data{4,ind} = t;
                data{5,ind} = v_ode;
                this.setOdeSave(ind +1,data);
            end
            
        end
        
        function clockEdgeTimes = findClockEdges(this,t)
            
            clkEdgeOptions = this.tbOptions.clockEdgeSettings;
            
            %produce a clock trace and find where the clock edge crosses
            %between 'low' percentage and 'high' percentage of vdd.
            clk = this.clockSource.V(t./this.tbOptions.capScale);
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
                if(abs(edges(i)-edges(i+1)) < clkEdgeOptions.clkEdgeSpread*this.tbOptions.capScale)
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
        
        function [index] = findNextBracket(this,Vin,t,edgeTimes)
            intHigh = 0;
            intLow = 100;
            
            currentInteval = this.dinInterval;
            dataTransitionTimes = linspace(currentInteval(1),currentInteval(2),...
                this.tbOptions.numParallelCCTs);
            edgeIndex = find(t < edgeTimes);
            if(isempty(edgeIndex))
                vNode = this.nodeMask(end);
            else
                vNode = this.nodeMask(edgeIndex(1));
            end
            threshold = this.tbOptions.digitalSimOptions.threshold;
            stepOut   = this.tbOptions.stepOut;
            Vlow = this.tbOptions.vdd*threshold;
            Vhigh = this.tbOptions.vdd*(1-threshold);
            numStatesPerCCT = length(Vin)/this.tbOptions.numParallelCCTs;
            
            index = zeros(1,2);
            
            for i = 1+stepOut:(this.tbOptions.numParallelCCTs-stepOut)
                %means settled low
                if(Vin(numStatesPerCCT*(i-1)+vNode) < Vlow)
                    %this checks to make sure that the trajectory that is
                    %allegedly going to go low actually goes low
                    if(Vin(numStatesPerCCT*(i-stepOut-1) + vNode) < Vlow)
                        intLow = dataTransitionTimes(i-stepOut);
                        index(1) = i-stepOut;
                    end
                end
                %means settled high
                if(Vin(numStatesPerCCT*(i-1)+vNode) > Vhigh)
                    %this checks to make sure that the trajectory that is
                    %allegedly going to go high actually goes high
                    if(Vin(numStatesPerCCT*(i+stepOut-1) + vNode) > Vhigh)
                        intHigh = dataTransitionTimes(i+stepOut);
                        index(2) = i+stepOut;
                        break;
                    end
                end
            end
            
            if(intLow == 100)
                error('bisection failed, all traces settled high')
            end
            
            if(intHigh == 0)
                error('bisection failed, all traces settled low')
            end
            
            this.setInterval([intLow,intHigh]);
            
        end
        
        function [startTime,newStates,VL,VH,diffStates,timeInt] = findNextStartStatesAndTime(this,Vin,t,index)
            
            nTraces = this.tbOptions.numParallelCCTs;
            numStatesPerCCT = length(Vin(1,:))/nTraces;
            [nSlices,~] = size(Vin);
            maxDiffStates = zeros(numStatesPerCCT,nSlices);
            span = linspace(1,nTraces,nTraces);
            degree = this.tbOptions.polyFitDegree;
            error = this.tbOptions.polyFitError;
            
            for i = 1:nSlices
                slice = Vin(i,:);
                slice = reshape(slice,[numStatesPerCCT,nTraces]);
                maxDiffStates(:,i) = slice(:,index(1))- slice(:,index(2));
                A = zeros(numStatesPerCCT,1);
                
                for j = 1:numStatesPerCCT
                    A(j) = polyval(polyfit(span,slice(j,:),degree),t(i));
                end
                Acompare = repmat(A,[1,nTraces]);
                if(norm(slice-Acompare,'fro') > error)
                    startTime = t(i);
                    timeInt = t(1:i);
                    diffStates = maxDiffStates(:,1:i);
                    break;
                end
            end
            
            newStates = zeros(numStatesPerCCT,nTraces);
            alpha = linspace(0,1,nTraces);
            VL = slice(:,index(1));
            VH = slice(:,index(2));
            
            for i = 1:nTraces
                newStates(:,i) = (1-alpha(i))*slice(:,index(1)) + alpha(i)*slice(:,index(2));
            end
            newStates = reshape(newStates,[1 numStatesPerCCT*nTraces]);
            newStates = newStates';
        end
        
        function [value,isterminal,direction] = stopCriteria(this,t,V)
            capScale = this.tbOptions.capScale;
            Vin = this.vfull(t,V);
            iss = this.is_src;
            Vin = Vin(~iss,:);
            Vin = reshape(Vin,size(V));
            vDot = this.dV_ode(t,V);
            [value,isterminal,direction] = stopToDigitalSim(t/capScale,Vin,vDot,this.tbOptions.digitalSimOptions);
        end
        
        function bisectionPlotV(this,plotHandles,t,Vin)
            
            plotOptions = this.tbOptions.plotOptions;
            
            %if the plotOptions are set to false then plot the default nodes
            %which are the nodes specified by the user which are decision nodes
            %for nested bisection and the approximation to beta plot
            
            %Reserved, the last plotHandle in the set of plot handles is for
            %plotting the approximation to beta
            
            [~,numStates] = size(Vin);
            statesPerCCT = numStates/this.tbOptions.numParallelCCTs;
            
            
            if(~plotOptions)
                
                nodes = this.nodeMask;
                
                for i = 1:length(plotHandles)
                    handles = plotHandles{i};
                    for j = 1:length(handles)
                        plot(handles(j),t,Vin(:,nodes(j)),'b--','linewidth',1.1)
                        for numCCTs = 2:this.tbOptions.numParallelCCTs-1
                            plot(handles(j),t,Vin(:,nodes(j)+(numCCTs-1)*statesPerCCT))
                        end
                        plot(handles(j),t,Vin(:,nodes(j)+(this.tbOptions.numParallelCCTs-1)*statesPerCCT),'r--','linewidth',1.1)
                    end
                end
            else
                
                numFigures = plotOptions.numFigures;
                nodesToPlotPerSubplot = plotOptions.nodesToPlot;
                %this flexibility allows the user to plot a specified number of
                %figures, where each figure has a number of subplots and for
                %each subplot the user can select what nodes to print. Hard
                %coded here is that this node is plotted for all parallel
                %circuits being simulated.
                
                for i = 1:numFigures
                    handles = plotHandles{i};
                    for j = 1:length(handles)
                        nodes =  nodesToPlotPerSubplot{j};
                        for k = 1:length(nodes)
                            for numCCTs = 1:this.tbOptions.numParallelCCTs
                                plot(handles(j),t,Vin(:,nodes(k)+(numCCTs-1)*statesPerCCT))
                            end
                        end
                    end
                end
                
            end
            drawnow
        end
        
        function bisectionPlotBeta(this,plotHandles,tbeta,beta)
            plotOptions = this.tbOptions.plotOptions;
            if(~plotOptions)
                plot(plotHandles,tbeta,log10(beta),'b')
            else
                %options you can think of would go here
            end
            drawnow
        end
        
        function plotHandles = bisectionPlotSetup(this,tspan,tCrit,tOff)
            plotOptions = this.tbOptions.plotOptions;
            
            %if there are no plot options create a subplot with each subplot
            %containing the node of interest for determining metastability
            %resolution and the clock source overlayed
            if(~plotOptions)
                nodes = this.nodeMask;
                
                plotHandles = cell(1,2);
                
                res = round((tspan(2)-tspan(1))/1e-12);
                t = linspace(tspan(1),tspan(2),res);
                clkSig = this.clockSource.V(t);
                tStart = this.tbOptions.transSettle;
                
                figure
                for i = 1:length(nodes)
                    outputPlots(i) = subplot(length(nodes),1,i);
                    plotTitle = sprintf('Node %d trace',nodes(i)+length(this.sources));
                    title(plotTitle);
                    ylabel('Output [V]');
                    set(gca,'fontsize',18);
                    axis([tStart-tOff,tspan(2)-tOff,0,1]);
                end
                xlabel('Time [s]') %make the last subplot contain the time xlabel
                %configure the subplots to add what is plotted on it
                set(outputPlots,'Nextplot','add');
                
                for i = 1:length(outputPlots)
                    plot(outputPlots(i),0,0,'linewidth',1.2) %ghost point for legend purposes
                    plot(outputPlots(i),t-tOff,clkSig,'-g','linewidth',3);
                    plot(outputPlots(i),ones(1,100)*tCrit-tOff,linspace(0,1),'m--','linewidth',2)
                end
                
                figure
                betaPlot = subplot(1,1,1);
                title('log_{10}||\beta(t)||_2 vs time');
                xlabel('Time [s]');
                ylabel('log Gain [V/s]');
                set(gca,'fontsize',18);
                set(betaPlot,'Nextplot','add');
                plotHandles{1} = outputPlots;
                plotHandles{2} = betaPlot;
            else
                
                %TODO fill in code to allow user configurable plot options
                % something along the lines of how many figures, how many nodes
                % to plot per subplot and finally an array of those nodes to plot per subplot
                % also add user configurability for linewidth and colors to
                % plot
                
                %             numFigures = plotOptions.numFigures;
                %             numSubplots = plotOptions.numSubplots;
                %             numNodesPerSubplot = plotOptions.nodesPerSubplot;
                %             nodesToPlotPerSubplot = plotOptions.nodesToPlot;
                %             plotTitles = plotOptions.titles;
                
            end
            drawnow
        end
        
        % compute the voltage for all circuit nodes
        % v_ode: each column for a time
        function vfull = vfull(this,t,v_ode)
            assert(length(t)==size(v_ode,2));
            
            v_odeShape = size(v_ode);
            
            v_ode = reshape(v_ode,[],this.tbOptions.numParallelCCTs);
            % voltages of sources at time t
            vs = this.eval_vs(t./this.tbOptions.capScale);
            % assign voltage source
            iss = this.is_src;
            vfull = zeros(this.circuit.nodeNum,v_odeShape(2)); %the full vector for one circuit
            vfull(this.map,:) = vs; % fill in voltage sources map to circuit port
            vfull = repmat(vfull,this.tbOptions.numParallelCCTs,1);   %get num parallel cct copies
            vfull = reshape(vfull,[],this.tbOptions.numParallelCCTs);
            vfull(~iss,:) = v_ode;
            
            %assign data input voltages
            interval = this.dinInterval;
            settleLow = interval(1);
            settleHigh  = interval(2);
            %get the array of delays
            delays = linspace(settleLow,settleHigh,this.tbOptions.numParallelCCTs);
            din = this.inputSource.V(t./this.tbOptions.capScale,delays./this.tbOptions.capScale);
            %place the differeing input voltages into vfull
            vfull(this.inputSourceIndex,:) = din;
            
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
                        vfull(vSrcMap(1),:) = max(min(v(2,:) + src.V(t./this.tbOptions.capScale),this.tbOptions.vdd),0);
                    end
                end
            end
            
            
            %if AD is true, then make vfull AD variables
            if(this.tbOptions.isAD)
                vfull = gradientinit(vfull);
            end
        end
        
        
        function setInterval(this,newInterval)
            this.dinInterval = newInterval;
        end
        
    end
end
