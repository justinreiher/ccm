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

classdef testbench < handle
  properties (GetAccess='public', SetAccess='private');
    % a circuit to be simulated
    circuit=[];
    % voltage sources
    sources={};
    % mapping between circuit ports and voltage sources
    map=[]; 
  end

  methods
    function this = testbench(circuit,ports,sources,do_opt)
      if(nargin<1), error('not enough parameters'); end
      if(nargin<2), ports = {}; end
      if(nargin<3), sources = {}; end
      if(nargin<4||isempty(do_opt)),do_opt=true;end
      if(length(ports)~=length(sources))
        error('incorrect size of sources and ports');
      end
      c = circuit;
      if(do_opt), c = c.flatten.vectorize; end % do optimization for performance
      this.circuit = c; 
      this.circuit = circuit;
      this.sources = sources;
      this.map = circuit.find_port_index(ports);
      if(any(this.map==0))
        error('can not find ports'); 
      end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Functions for simulation
    
    % This function simulate the circuit
    function [t,v,dv] = simulate(this,tspan,v0,opts)
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
        [t, v_ode] = integrator(@(t, v_ode)(this.dV_ode(t, v_ode)), tspan, v0_ode, opts); 
      end
      
      % set the value of inputs
      v = this.vfull(t,v_ode')'; % each row for a time point
      if(nargout>2) 
        dv_ode = this.dV_ode(t,v_ode')'; 
        tt = t; tt(end+1) = 2*t(end) - t(end-1);
        vIn = v(:,~isOde); 
        vIn(end+1,this.map) = this.eval_vs(tt(end))';
        dv_vs = (vIn(2:end,:) - vIn(1:end-1,:))./repmat(tt(2:end)-tt(1:end-1),1,size(vIn,2)); 
        dv = zeros(size(v));
        dv(:,isOde) = dv_ode;
        dv(:,~isOde) = dv_vs;
      end
    end % simulate 
    
    % This function is a wrapper of circuit.dV by setting input values 
    % as source.V(t) for the Matlab integrator.
    function dv = dV_ode(this,t,v_ode) 
      vfull = this.vfull(t,v_ode);
      dv = this.circuit.dV(vfull); % each column for a time point
      iss = this.is_src;
      dv = dv(~iss,:); 
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
  end
end
