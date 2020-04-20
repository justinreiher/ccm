% testbench properties
%   circuit: the circuit to be tested
%   sources: voltage sources for each circuit inputs
%   map:    mapping between voltage sources and circuit inputs
%   modelDictionary: The model dictionary object
%   defaultModel:    The default model used in this testbench
%   tbOptions:       The options used for this testbench
%   odeSave:         Cell array, populated if testbench is set in debug
%                    mode via the testbench options
%   matrixMap:       A matrix map which distributes capacitance and
%                    currents to the appropriate nodes
%   simVoltage:      The saved voltage vector after a simulation run, used
%                    for the provided plotting functions
%   simTime:         The saved simulation time points, corresponds to the
%                    simVoltage points.
%
% constructor
%   testbench(circuit,ports,sources,defaultModel,...
%             defaultModelProcess,testbenchOptions); 
%     circuit: an object of "circuit" class.
%     ports:    an cell of "node" objects
%     sources: an cell of "circuit" objects which must implement the V(t) function. 
%          the size of ports and sources must be the same. 
%          the voltage source is applied to the corresponding port
%     
%     optional inputs:
%       defaultModel: The global default model to use for simulation
%       defaultModelProcess: The global transistor process used in
%       simulation, i.e. PTM 180nm,PTM 45nmHP, ect
%       testbenchOptions: Options for the simulator, can be provided as a
%       set of key value pair of settings or bundled as a struct
%    
% methods
%   It provides functions for simulation
%   [t,v,debug] = simulate(tspan,v0,opts): this function simulate the circuit. 
%     tspan:  the range of simulate time 
%     v0:     the initial circuit state  
%     opts:   opt.integrator is a function handle for ODE integrator.
%             it is also used in the integrator function.
%     t:      simulated time point
%     v:      voltage of circuit nodes, each row for a time 
%     debug:  Capacitance, voltage, time and derivatives of node voltages,
%             each row for a time stored in a cell. Only provided if debug
%             flag is set to true in tbOptions.
%
%   dv_ode = dV_ode(t,v_ode): construct ODE for ODE nodes
%     t:      used to evaluates voltage of sources
%     v_ode:  voltage of ODE nodes
%     dv_ode: time derivative of ODE voltages
%
%   tbOptions = packageTestbenchOptions: Method to package testbench options into a
%   struct
%
%   tbOptions = getTestbenchOptions: Method to return current testbench options
%  
%   setTestbenchOptions:  Method to apply new set of testbench options to
%                         this testbench
%
%   setSimTimeVoltage:  Method to apply user specified time and voltage to
%                       the testbench. Allows the user to supply a
%                       previously saved trajectory to the testbench.
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
%   n  = dim: the number of dimensions of the system dyanamics
%
%   ind = find_port_index: wrapper of circuit's "find_port_index" funcition.
%
%   n = inputNum: the number of input sources
%
%   is = inputs: return the voltage sources (same as sources)
%
%   c = dut: return the circuit
%
%   ind = inputNodes: returns the map of signals
%
%   currents = getIdsDevices(query,time,v): returns the current for the set
%   of devices specified by the query. Requires specifying all the nodes
%   associated with the device. More than one device can be queried at a
%   time. The time and v parameters are optional, if they are not
%   provided, they are retrieved from the saved simulation data.
%
%   explainNodes(nodes,time, vSupplied): returns a 2x1 plot for each node
%   in nodes. Each plot show the contribution of all devices providing
%   current into the node in question. The voltage on the node in question
%   is also plotted. The second plot shows the sum contribution of all the
%   currents aswell as the voltage on the node in question
%
%   setOdeSave: a function which updates the odesSave structure when the
%   testbench is in debug mode
%
%   odeSave = getOdeSave: returns the object where all the ode data is
%   saved
%
%   iss = is_src: returns a boolean vector where a 1 is placed if the node
%   is a source node.
%
%   maps = computeMapMatrix: returns a cell array which contains the
%   mapping information for each unique kind of circuit element present in
%   the circuit
%

% Justin Reiher 2020: Change log
%   Added a model dictionary object which holds information on which models
%   are available, such as: MVS, EKV, Interpolated
%   Additional optional options allow a user to specify the model to use
%   and the testbench options.
%
% The testbench options are stored as a struct with the following key,value pairs:
%   capModel: full, means a capacitance model where miller capacitances are
%   included. Future work could include wire capacitance
%               gnd, means all capacitances are referenced to ground
%   capScale: This is a parameter to allow scaling of time,used to help
%   numerical intergation
%   vdd: This is the global voltage supply applied to non-specified vdd
%   nodes. Legacy from some previously defined circuit models
%   temp: The temperature the circuit is to be simulated at. Only
%   applicable if the model takes into account temperature
%   numParallelCCTs: Allows the user to specify how many copies of the
%   circuits are being simulated. If more than one, then different initial
%   conditions can be applied to the different copies of the circuit
%   debug: This is a boolean variable which allows if set to true saves all
%   the data as the numerical integration takes place
%
%   The testbench options can also be specified as key,value pairs. This is
%   not recommended, but a function exsists to take the list of key,value
%   pairs and package them into a struct if constructed in this way.
%
% It is recommended to add new testbench options in the Matlab convection
% of a 'string', value pair. 

classdef testbench < handle
  properties (GetAccess='public', SetAccess='private');
    % a circuit to be simulated
    circuit=[];
    % voltage sources
    sources={};
    % mapping between circuit ports and voltage sources
    map=[];
    %the model dictionary that knows about all the models to compute dv for
    %the various device types know for the simulator.
    modelDictionary = [];
    %sets the default model.
    defaultModel = [];
    %sets simulation options, i.e. model accuracy
    tbOptions = [];
    %debug info if set
    odeSave = [];

    %Saved simVoltage and simTime
    simVoltage = [];
    simTime    = [];
    
  end
  properties (GetAccess = 'protected', SetAccess = 'private')
      circuitMap;
      %total circuit mapping
      matrixMap = [];
  end
  
  methods
      % to create a testbench you need to provide a circuit description,
      % the associated ports and a model dictionary to load the appropriate
      % default model to use for the simulation. The default tbOptions are
      % such that it works for the original default 180nm process which ran
      % at 1.8V. The other default is that all capacitances go to ground.
    function this = testbench(circuit,ports,sources,...
            defaultModel,defaultModelProcess,varargin)
       tbOptions = [];
      if(nargin<1), error('not enough parameters'); end
      if(nargin<2), ports = {}; end
      if(nargin<3), sources = {}; end
      if(nargin<4), defaultModel = 'interp'; defaultModelProcess = 'PTM 180nm'; end
      if(nargin<5), tbOptions = struct('capModel','gnd','capScale',1,'vdd',1.8,'temp',298,'numParallelCCTs',1,'debug',false); end
      if(length(ports)~=length(sources))
        error('incorrect size of sources and ports');
      end

      if(isempty(tbOptions))
          tbOptions = this.packageTestbenchOptions(varargin);
      end
      this.modelDictionary = modelDictionary;
      this.defaultModel = modelDictionary.getModel(defaultModel,defaultModelProcess,tbOptions);
      tbOptions.numNodes = circuit.nodeNum;
      tbOptions.numStates = circuit.nodeNum - length(sources);
      this.tbOptions = tbOptions;
      c = circuit;
      c = c.flatten.vectorize; % jr: circuit must be flattened and vectorized
      %if the defaultModel class cannot compute vectorized quantities it
      %needs to evaluate them one at time. Not the responsibility of the
      %testbench.
      this.circuit = c; 
    %  this.circuit = circuit;
      this.sources = sources;
      this.map = circuit.find_port_index(ports);
      this.modelDictionary = modelDictionary();
      this.circuitMap = this.circuit.getCircuitMap;
      if(any(this.map==0))
        error('can not find ports'); 
      end
       
      % computes the matrix map to simulate the number of specified
      % parallel circuits
      maps = this.computeMapMatrix;
      M = cell(1,length(maps));
      for i = 1:length(maps)
        M{i} = sparse(kron(eye(this.tbOptions.numParallelCCTs),maps{i}));
      end
      this.matrixMap = M;
       
       if(this.tbOptions.debug)
           odeSave = cell(1,5000);
           odeSave{1} = 1;
           this.odeSave = odeSave;
       end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Functions for simulation
    
    % This function simulate the circuit
    function [t,v,debug] = simulate(this,tspan,v0,opts)
      n = this.circuit.nodeNum;
      % initial values
      if(nargin < 3) 
        v0 = zeros(n,1); % initialize all nodes to 0 
      elseif(iscell(v0)) % user specified a few nodes to initialize 
        v00 = v0;  
        v0 = zeros(n,1); % set the others to 0 
        v0(v00{1}) = v00{2}; 
      end 
      % Integrator
      if(nargin < 4), opts = odeset(); end; 
      if(isfield('Integrator', opts)) 
        integrator = opts.Integrator; 
      else 
        integrator = @ode45; 
      end 

      % simulate ODE nodes
      isOde = ~this.is_src;
      v0_ode = v0(isOde);  
      if(isempty(v0_ode))
        t = reshape(tspan, [], 1);
        v_ode = zeros(length(t), 0);
      else
        tspan = this.tbOptions.capScale*tspan;
        [t, v_ode] = integrator(@(t, v_ode)(this.dV_ode(t, v_ode)), tspan, v0_ode, opts); 
      end
      
      % set the value of inputs
      v = this.vfull(t/this.tbOptions.capScale,v_ode')'; % each row for a time point
      t = t/this.tbOptions.capScale;
      this.simVoltage = v;
      this.simTime = t;
      if((nargout > 2)&&this.tbOptions.debug)
          debug = this.odeSave;
      end
    end % simulate 
    
    % This function is a wrapper of circuit.dV by setting input values 
    % as source.V(t) for the Matlab integrator.
    function dv = dV_ode(this,t,v_ode)
      
      vodeShape = size(v_ode);
      vfull = this.vfull(t,v_ode);
      circuitMapping = this.circuitMap;

      %%%%%%
      %The circuit map returns a cell array where each entry contains is of
      %the form {ElementProperties,Map}, the elements are a cell array of all the
      %same kind of elements and the Map is the map that corresponds to the
      %appropriate mapping for those elements.
      %%%%%
      [~,numUniqueElements] = size(circuitMapping);
      Itot = cell(1,numUniqueElements);
      Ctot = cell(1,numUniqueElements);
      Ms   = this.matrixMap;
      for i = 1:numUniqueElements
          %get the circuit element grouping
          elementGroup = circuitMapping{i};
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
      I = I(~isSource);
      C =  this.tbOptions.capScale * C;
      
      if(strcmp(this.tbOptions.capModel,'gnd'))
          C = C(~isSource);
          dv = I./C;
      else

          %C here is the diagonal matrices of the sparse matrix I want to have,
          %MATLAB doesn't offer a convienient way to create a block
          %diagonal sparse matrix other than by the following:
          B = mat2cell(sparse(C),numCircuitNodes*ones(1,this.tbOptions.numParallelCCTs),numCircuitNodes);
          C = blkdiag(B{:});
          C = C(~isSource,~isSource);
          dv = C\I;
      end

      dv = reshape(dv,vodeShape);
      % It the debug flag is set, then save all ode data at each time point,
      % which includes:
      %     C  - The computed capacitance matrix at time t
      %     I  - The computed current at time t
      %     dv - The computed derivative at time t
      %     t  - Save the time point t
      % This debug feature is not memory efficient, and stores the above 4
      % quantities in a 4,1 cell which is stored in a cell that grows at
      % each time point. Useful to observe what is happening for debugging
      % purposes but likely to dramatically slow down simulation time
      if(this.tbOptions.debug)
          odeSave = this.odeSave;
          ind = odeSave{1};
          data = cell(4,1);
          data{1} = C;
          data{2} = I;
          data{3} = dv;
          data{4} = t;
          odeSave{ind+1} = data;
          odeSave{1} = ind + 1;
          this.odeSave = odeSave;
      end
    end
    
    %%%
    %Function to package testbench options into a struct if provided as a
    %key,value pairing. If the options are already provided as a struct,
    %return the struct and if there is an odd number of parameters produce
    %an error message notifying the user of the mismatch.
    %%%
    function options = packageTestbenchOptions(this,tbOptions)
        if(isstruct(tbOptions{1}))
            options = tbOptions{1};
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
        end
            
    end
    
    %setter and getter methods
    function tbOptions = getTestbenchOptions(this)
        tbOptions = this.tbOptions;
    end
    
    function setTestbenchOptions(this,tbOptions)
        this.tbOptions = tbOptions;
    end
    
    function setSimTimeVoltage(this,time,voltage)
        this.simTime = time;
        this.simVoltage = voltage;
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Functions for coho

    % This function computes LDI models for input signals specified by Brockett annuli
    % region: current reachable region (an object of shape class) 
    % states: the current region (1-4) for Brockett annulus of each signal.
    % vars:   index of interested nodes
    function [c,err] = dV_ldi(this,region,states,vars,varargin)
      if(nargin<4||isempty(vars)), vars = 1:this.circuit.nodeNum; end
      vars = reshape(vars,1,[]);
      [c,err] = this.circuit.dV_ldi(region,vars,varargin{:});
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
      n = this.circuit.nodeNum;
    end
    % The name of the circuit
    function str = name(this)
      str = class(this.circuit);
    end
    % design under test (circuit)
    function c = dut(this)
      c = this.circuit;
    end
    % The number of input signals 
    function n = inputNum(this)
      n = length(this.sources);
    end
    % voltage source of input signals
    function i = inputs(this,ind)
      if(nargin<2||isempty(ind))
        i = this.soruces;
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
    
    % Function which returns the current flowing through a device or pair
    % of nodes determined by the parameter query
    
    function currents = getIdsDevices(this,query,time,v)
        
        cct = this.circuit;
        cctMap = cct.maps;
        
        if(nargin < 3) v = this.simVoltage'; end
        
        if(length(v(:,1)) ~= cct.nodeNum)
            v = this.vfull(time,v);
        end
            
        [numTerminals,numDev] = size(query);
        
        nIdentifier  = logical(zeros(1,numDev));
        currents      = zeros(numDev,length(v));
        wid           = cell(1,numDev);
        rlen          = cell(1,numDev);
        model         = cell(1,numDev);
        devMap        = zeros(4,numDev);
        
        Vdev = zeros(numTerminals,length(v),numDev);
        
        for j = 1:numDev
            for i = 1:length(cctMap)
                if(cctMap{i} == query(:,j))
                    device   = cct.getDevicesByIndex(i);
                    wid{j}   = device{1}.params.wid;
                    rlen{j}  = device{1}.params.rlen;
                    model{j} = device{1}.params.model;
                    devMap(:,j) = cctMap{i};
                    if(strcmp(device{1}.subinfo.device,'nmos'))
                        nIdentifier(j) = 1;
                    end
                    Vdev(:,:,j) = v(cctMap{i},:);
                end
            end
        end
        
        vectorizedParamsN = struct('wid',{wid(nIdentifier)},'rlen',{rlen(nIdentifier)},'model',{model(nIdentifier)},'deviceMap',{devMap(:,nIdentifier)});
        vectorizedParamsP = struct('wid',{wid(~nIdentifier)},'rlen',{rlen(~nIdentifier)},'model',{model(~nIdentifier)},'deviceMap',{devMap(:,~nIdentifier)});
        
        Vn = Vdev(:,:,nIdentifier);
        Vp = Vdev(:,:,~nIdentifier);
        
        if(~isempty(Vn))
            Idn = this.defaultModel.I('nmos',vectorizedParamsN,...
                                          Vn,this.tbOptions);
        end
        
        if(~isempty(Vp))
            Idp = this.defaultModel.I('pmos',vectorizedParamsP,...
                                          Vp,this.tbOptions);
        end
        
        nCounter = 1;
        pCounter = 1;
        
        for i = 1:numDev
            if(nIdentifier(i))
                currents(i,:) = Idn(1+(nCounter-1)*4,:);
                nCounter = nCounter +1;
            else
                currents(i,:) = Idp(1+(pCounter-1)*4,:);
                pCounter = pCounter + 1;
            end
        end
        
    end
    
    % Function which "explains" the nodes in question of the circuit over
    % the simulation time span [0,T], for every node the function returns a:
    %       Plot which shows the individual components of the
    %       currents flowing into the node(s) in question by any device
    %       connected to the node in question along with the voltage
    %       present on that node
    %
    %       Plot which sums the current from all the individual
    %       contributions from the node(s) along with the voltage present
    %       on that/those node(s).
    %
    % The user can supply their own simulation time and voltage vector or
    % after a simulation is complete the explainNodes function takes the
    % simTime and simVoltage from the most recently run simulation. If a
    % simulation has not been performed an error is returned.
    
    function explainNodes(this,nodes,time,vSupplied)
        cct = this.circuit;
        if(nargin < 3) 
            vSupplied = this.simVoltage';
            time = this.simTime;
            if(isempty(vSupplied) || isempty(time))
                error('Testbench has no simulation data to display, consider running a simulation');
            end
        end
        
        [numNodes,numTimePoints] = size(vSupplied);
        if(numNodes < cct.nodeNum)
            vSupplied = this.vfull(time,vSupplied);
        end
        
        for i = 1:length(nodes)
            figure
            subplot(1,2,1)
            hold on
            [devMaps,names] = cct.getDevicesConnectedToNode(nodes(i));
            currents = this.getIdsDevices(devMaps,time,vSupplied);
            plot(time,currents,'-','linewidth',1.2)
            ylabel('Current [A]')
            yyaxis right
            plot(time,vSupplied(nodes(i),:),'m-','linewidth',1.2)
            ylabel('Voltage [V]')
            names{end+1} = strcat('Voltage Node - ',cct.find_node_name_short(nodes(i)));
            xlabel('Time [s]')
            
            legend(names)
            set(gca,'fontsize',20)
            
            subplot(1,2,2)
            hold on
            plot(time,sum(currents),'-','linewidth',1.2)
            ylabel('Current [A]')
            yyaxis right
            plot(time,vSupplied(nodes(i),:),'m-','linewidth',1.2)
            ylabel('Voltage [V]')
            xlabel('Time [s]')
            
            legend(strcat('Net current Node - ',cct.find_node_name_short(nodes(i))),strcat('Voltage Node - ',cct.find_node_name_short(nodes(i))));
            set(gca,'fontsize',20)
        end
    
    end

    function setOdeSave(this,ind,data)
        odeSaveCopy = this.odeSave;
        odeSaveCopy{1} = ind;
        odeSaveCopy{2} = data;
        this.odeSave = odeSaveCopy;
    end
    
    function odeSave = getOdeSave(this)
        odeSave = this.odeSave;
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
      n = this.circuit.nodeNum;
      iss = false(n,1); 
      iss(this.map) = true;
    end
    
    %Function which computes the circuit matrix maps for the device types
    %in the circuit (nmos,pmos,vsrc, ...ect). Computed once and stored to
    %be used to compute the current Icct = Mdev*Idev

    function maps = computeMapMatrix(this)
        
        devMaps = this.circuitMap;
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
  end
    
    methods(Access='private')
        % compute the voltage for all circuit nodes
        % v_ode: each column for a time
        
        function vfull = vfull(this,t,v_ode)
            assert(length(t)==size(v_ode,2));
            % voltages of sources at time t
            vs = this.eval_vs(t);
            % assign voltage source
            iss = this.is_src;
            vfull = zeros(this.circuit.nodeNum,size(v_ode,2));
            vfull(~iss,:) = v_ode;
            vfull(this.map,:) = vs; % map to circuit port
            
            %Also need to add any internal voltage sources that may exist
            
            circuitMapping = this.circuitMap;
            %%%%%%
            %The circuit map returns a cell array where each entry contains is of
            %the form {ElementProperties,Map}, the elements are a cell array of all the
            %same kind of elements and the Map is the map that corresponds to the
            %appropriate mapping for those elements.
            %%%%%
            [~,numUniqueElements] = size(circuitMapping);
            
            for i = 1:numUniqueElements
                dev = circuitMapping{i};
                devParams = dev{1};
                devMap = dev{2};
                if(strcmp(devParams.deviceType, 'vsrc'))
                    internalVsrc = dev{3};
                    for j = 1:length(internalVsrc)
                        src = internalVsrc{j};
                        vSrcMap = devMap(:,j);
                        v = vfull(vSrcMap,:);
                        vfull(vSrcMap(1),:) = v(2,:) + src.V(t);
                    end
                end
            end
            
        end
  end
end
