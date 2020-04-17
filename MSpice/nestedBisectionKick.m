% nestedBisection is a superclass of nestedBisectionKick, has no properties
%
% constructor
%   nestedBisectionKick(synchronizerCCT,ports,sources,dinInterval,dinSource,...
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
% methods
%   It provides functions for simulation
%   [t,v,dv] = simulate(filename,tspan,tCrit,kickData,v0,opts): this function simulate the circuit. 
%     filename: a string where all the nested bisection data is stored in a
%     .mat file
%     tspan:  the range of simulate time 
%     tCrit:  the user specified critical time at which the synchronizer
%     must be at a settled value, used as a criteria to stop the simulation
%     kickData: a struct of the format:
%                   kickTime - the time in the simulation at which the kick
%                   occurs
%                   kickPercentage - the percentage amount by which the
%                   kicked node/nodes of interest are "kicked"
%                   kickNodes - the nodes that are modified at time kickTime
%                   by kickPercentage
%     v0:     the initial circuit state  
%     opts:   opt.integrator is a function handle for ODE integrator.
%             it is also used in the integrator function.
%     t:      simulated time point
%     v:      voltage of circuit nodes, each row for a time 
%     dv:     derivatives of node voltages, each row for a time  
%


classdef nestedBisectionKick < nestedBisection
  
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
      
    function this = nestedBisectionKick(synchronizerCCT,ports,sources,dinInterval,dinSource,...
            clockSource,nodeMask,modelDictionary,defaultModel,defaultModelProcess,varargin)
      
      this = this@nestedBisection(synchronizerCCT,ports,sources,dinInterval,dinSource,...
            clockSource,nodeMask,modelDictionary,defaultModel,defaultModelProcess,varargin);
            
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Functions for simulation
    
    % This function simulate the circuit
    function [bisectionRuns] = simulate(this,filename,tspan,tCrit,kickData,v0,opts)
        
        if(ischar(filename)||isstring(filename))
            fprintf('starting nested bisection and saving data to filename: %s ...\n',filename);
        else
            error('filename needs to be a string or array of characters')
        end
        if(isempty(kickData))
            error('Kick data needs to be specified as struct: kickTime,kickPercentage,kickNodes')
        end
        
        kickTime = kickData.kickTime*this.tbOptions.capScale;
        kickPercentage = kickData.kickPercentage;
        kickNodes = kickData.kickNodes;
        
        n = this.circuit.nodeNum;
        tspan = tspan*this.tbOptions.capScale;
        tCrit = tCrit*this.tbOptions.capScale;
        % initial values
        if(nargin < 6)
            v0 =zeros(n,1); % initialize all nodes to zero
        elseif(iscell(v0)) % user specified a few nodes to initialize
            v00 = v0;
            v0 = zeros(n,1); % set the others to 0
            v0(v00{1}) = v00{2};
        end
        % Integrator
        %set the default integrator options if none are specified
        if(nargin < 7), opts = this.tbOptions.integratorOptions; end
        if(isfield('Integrator', opts))
            integrator = opts.Integrator;
        else
            integrator = this.tbOptions.integrator;
        end
        
        % simulate ODE nodes
        isOde = ~this.is_src;
        v0_ode = v0(isOde);
        v0_ode = repmat(v0_ode,[this.tbOptions.numParallelCCTs,1]); %set up n parallel circuits to run concurrently
        
        initialInterval = this.dinInterval;
        
        if(isempty(v0_ode))
            t = reshape(tspan, [], 1);
            v_ode = zeros(length(t), 0);
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

            fprintf('Starting nested bisection...\n')
            tic()
            while(~stop)
                
                if(this.tbOptions.debug)
                    ind = 1;
                    data = cell(5,5000);
                    this.setOdeSave(ind,data);
                end
                %if the start time is less than the time at which the kick
                %occurs, then we simulate upto the kickTime or less if the
                %synchronizer settles
                if(startTime < kickTime)
                    [tBeforeKick, vBeforeKick,timeEvent,vEvent,indexEvent] = integrator(@(tBeforeKick, vBeforeKick)(this.dV_ode(tBeforeKick, vBeforeKick)), [startTime,kickTime], v0_ode, opts);
                    %Here we check to see if we actually simulated to the
                    %time of the kick, if we have then we apply the kick,
                    %if not then we resume bisection as normal. If we have
                    %simulated upto the time of the kick, we apply the kick
                    %and keep simulating.
                    v_ode = vBeforeKick;
                    t = tBeforeKick;
                    if(tBeforeKick(end) >= kickTime)
                        icCont = reshape(vBeforeKick(end,:),[numStatesPerCCT, this.tbOptions.numParallelCCTs]);
                        icCont(kickNodes,:) = kickPercentage*icCont(kickNodes,:);
                        icCont = reshape(icCont,[1, numStatesPerCCT*this.tbOptions.numParallelCCTs]);
                        [tAfterKick,vAfterKick] = integrator(@(tAfterKick,vAfterKick)(this.dV_ode(tAfterKick,vAfterKick)),[kickTime,stopTime], icCont,opts);
                        v_ode = [vBeforeKick;vAfterKick];
                        t = [tBeforeKick;tAfterKick];
                    end
                    %If we have progressed beyond the time of the kick,
                    %then we continue simulating as normal.
                else
                    [t,v_ode,timeEvent,vEvent,indexEvent] = integrator(@(t,v_ode)(this.dV_ode(t,v_ode)),[startTime,stopTime],v0_ode,opts);
                end
                
                indAfterSettle = t > this.tbOptions.transSettle*this.tbOptions.capScale;
                t = t(indAfterSettle);
                tEnd = t(end);
                v_ode = v_ode(indAfterSettle,:);
                this.bisectionPlotV(plotHandles(1),t./this.tbOptions.capScale - tOff,v_ode)  
                
                
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
                vTrajH = v_ode(:,1+numStatesPerCCT*(index(2)-2):numStatesPerCCT*(index(2)-1))';    
                vTrajL = v_ode(:,1+numStatesPerCCT*(index(1)):numStatesPerCCT*(index(1)+1))';
                
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
   
  end
end
