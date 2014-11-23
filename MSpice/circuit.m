%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class defines an abstract circuit. It provides template functions based
% on the concept that a circuit is a connected set of sub-circuits/elements. 
% It represents a circuit as a hierarchy tree, traverse the tree, call nodes' 
% function, merge results for the circuit. Leaf nodes are called as leaf-circuits, 
% which must implement and override the template functions by 'circuit' class.
% With leaf-circuits provided, a complicated circuit can be construct easily 
% will all function ready from this 'circuit' superclass. 
%
%------------------------------------------------------------------------------- 
% For using the circuit for simulation/verification/small-signal-analysis, etc.
%------------------------------------------------------------------------------- 
% The circuit class provides the following methods and members automatically:
%   circuit properties:  
%     name                                 % name of the circuit
%     nodeNum:                             % total number of nodes
%     ports:                               % visible(may not be all) nodes 
%     finalized;                           % status 
%     flattened;                           % status
%     vectorized;                          % status
%   constructors:  N/A 
%   methods: 
%     I. Simulation Interface
%       1. i = I(V)                        % current out of each node  
%       2. c = C(V)                        % capacitance
%       3. dv = dV(V)                      % ODE model (dv/dt = i/c)
%       4. v = V(t)                        % for voltage source
%     II. Verification Interface
%       1. [c,err] = I_ldi(region,vars)    % linear differential inclusion (LDI) model of currents
%          region:  interested region, objects of "Shape" class. 
%          vars:    index of interested nodes, all nodes by default 
%          i \IN c'*[V;1]+/-err
%       2. [c,err] = dV_ldi(region,vars)   % LDI models of dv/dt 
%          dv \IN c'*[V;1]+/-err
%          NOTE: currently we assume capacitance is independent of voltage for LDI models. 
%     III. Small Signal Analysis Interface
%       1. j = Jac(V)                      % Jacobian matrix
%       2. di = dIdV(V)                    % di/dv, for Jacobian matrix
%     IV. Query Functions:  a circuit may not support all interfaces above, which can be queried by
%       1. is = ifc_simu                   % Support interface I. 
%          a. is = is_vsrc                 % support I.4 if is voltage source; I.1/2/3 otherwise
%       2. is = ifc_verify                 % Support interface II.
%       3. is = ifc_ssan                   % Support interface III. 
%     V.  Optimization Functions: functions for improving simulation performance 
%       1. this = flatten                  % Flatten the circuit tree
%       2. this = vectorize                % Group the sub-circuits by class, used for I_vec/C_vec
%         NOTE: 'flatten' reduce the tree depth to be 1, and updates sub-circuits/mapping correspondingly 
%               'vectorize' groups sub-circuits, makes I call I_vec, replacing the default I_seq
%         NOTE: usually, c.flatten.vectorize gives the best performance for I/C for simulation. 
%         NOTE: these functions change the internal circuit tree, to save the original circuit,
%               please use copy(obj) before the function call. 
%     VI. Utility Functions 
%       1. is = is_leaf                    % Leaf circuit is the one without any sub-circuits
%       2. num = elemNum                   % Number of sub-circuits
%       3. ind = find_port_index(ports)    % Find the ID of a circuit port 
%       4. names = find_node_name(index)   % Return node names for given index
%       5. strs = print_circuit_tree(prefix)% Print the circuit tree structure   
%       6. strs = print_status(prefix)     % Print circuit and sub-circuit status
%       7. strs = print_nodes              % Description of all nodes
%
%   static methods
%     % Given a list of objects with same class, compute the current and capacitance. 
%     i = I_objs(Objs,V);  
%     c = C_objs(Objs,V);  
%   
%------------------------------------------------------------------------------- 
% For constructing complex circuits 
%------------------------------------------------------------------------------- 
% Users can define a subclass of 'circuit' for any complex circuits. 
% The subclass must provide (only) a public constructors, which follows
%   1. call superclass constructor, e.g. 
%      this = this@circuit('circuit_name'); 
%   2. compose the circuit by sub-classes using 
%      a. this.port = add_port(node('port_name'));  % add a new port to the circuit. 
%                                                   % port is an object of node class
%      b. e = add_element(element)                  % add an element(circuit) to the circuit. 
%                                                   % element is an object of circuit class
%      c. connect(this.port,e.port1,...)            % connect circuit ports and elements' ports
%   3. finalize the circuit
%      this.finalize;
%   NOTE: it's recommended to add circuit ports as public class members, such that ports can be easily 
%         accessed outside (e.g. used as elements for the connect command above). Otherwise, has to 
%         find in the 'ports' cell in the 'circuit' class.
%
%------------------------------------------------------------------------------- 
% For implementing leaf-circuits 
%------------------------------------------------------------------------------- 
% Leaf circuits can not reuse this class's template functions as there is no 
% sub-circuits in the leaf circuits. The leaf-circuit subclass must override the 
% following functions 
%   1. I     (if ifc_simu & ~is_vsrc)
%   2. C     (if ifc_simu & ~is_vsrc)
%   3. V     (if ifc_simu &  is_vsrc) 
%   4. I_ldi (if ifc_verify)
%   5. dIdV  (if ifc_ssan) 
%   6. Static.I_objs for vectorization (if for better performance)
%   7. Static.C_objs for vectorization (if for better performance)
%   8. vectorize_subtype_id (only if class have different nodeNum with configurations)
% All query functions return false by default, assuming leaf-circuits don't provide
% any of these functions. Query functions need to be enabled by subclass if provided. 
%
%------------------------------------------------------------------------------- 
% For understanding the implementation details
%------------------------------------------------------------------------------- 
% 1. Circuit Tree
%
%   We maintains the circuit tree by the private data 'elements' and 'maps'.
%   'elements' contains the list of sub-circuits, and 'maps' contains the mapping 
%   information between circuit nodes and sub-circuits nodes.
%
% 2. Flatten & Vectorize. 
%
%   'Flatten' re-constructs the circuit tree to make it flat, i.e. all elements are leaf-circuits. 
%   This helps to reduce the level of trees, thus the number of functional calls, and improve performance.  
%
%   'Vectorize' groups the elements by their class and merge all elements for vectorization. 
%   It uses 'grp_cids' to records elements IDs of the same type, and 'grp_maps' for vectorized mapping 
%   information for all such elements. These informations are used by functions 'I_vec/C_vec'.  
% 
%   Usually, c.flatten.vectorize provides the best performance for simulation. 
%   However, 'Vectorize' doesn't require prior 'Flatten', i.e. you can performance 'vectorize' without 
%   'flatten'. It's OK to do 'vectorize' before 'flatten' as well, but 'flatten' will destroy the internal 
%   data used by 'vectorize', thus need re-vectorize by users. 
%
%   When vectorized, 'I/C' function will call the vectorized function 'I_vec/C_vec'.
%   These function merges all function calls for elements of the same type to be one function calls
%   with vectorized parameters. Therefore, the performance is expected to be improved significantly. 
%   The one call is via the static methods 'I_objs/C_objs', which should be overridden by subclass.
%   Otherwise, the default template function from 'circuit' class is used, which calls elements' I/C
%   functions individually in a for loop, thus the same with I_seq, and no performance benefits. 
%
%   We don't support vectorized LDI model calculations, which would have little benefit of performance. 
%     
%------------------------------------------------------------------------------- 
% For further development 
%------------------------------------------------------------------------------- 
% MSpice is mainly for MOS based circuits. It has the following limitations
% 1. We use modified nodal analysis. But the voltage source is not supported by this 'circuit' class. 
%    The voltage source is treated as inputs of the circuit, and supported by the 'testbench' class. 
%    In other word, the DAE generated after MNA is converted to ODE by 'testbench' which are handled by 'circuit'. 
%    We may use 'supernode' or other technique to get ODE in this 'circuit' class directly.
% 2. We currently support first order ODE only, so inductor is not supported now. 
% 3. The ODE doesn't support time, i.e. vdot = f(v), independent of time t. That's the reason this class
%    doesn't not support time-dependent voltage/current sources. We may support it for simulation first, then verification.
% 

classdef circuit < matlab.mixin.Copyable %handle
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Class members
  %   Public: name, number of nodes, list of visiable ports
  %   Private: for circuit tree, vectorization, status, temp variables
  %
  properties (GetAccess='public', SetAccess='private');
    name='';         % circuit name
    nodeNum=0;       % # of circuit nodes
    ports={};        % visiable nodes
    %-- status --% 
    finalized=false;  % construction is complete 
    flattened=false;  % circuit tree flattened
    vectorized=false; % elements grouped for vectorization 
  end
  properties (GetAccess='protected', SetAccess='private');
    %-- circuit tree --%
    elements={};      % subcircuits 
    maps={};          % subcircuit nodes <-> circuit nodes
    conns={};         % connection info for computing 'maps' 
    %-- grouping for vectorization --% 
    grp_cids=[];      % class ID of each elements   
    grp_maps={};      % merged maps of all elements of the same class 
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Protected methods for composing larger circuits 
  %   constructor     % set name
  %   add_port        % set ports
  %   add_elements    % set elements
  %   connect         % set conns
  %   finalize        % set maps and nodeNumm
  methods(Access=protected)
    function this = circuit(name)
      if(nargin<1), name='anonymous'; end % TODO: remove?
      this.name = name;
    end

    function port = add_port(this,port)
      assert(~this.finalized); 
      for i=1:length(this.ports)
        if(this.ports{i}==port) % Is the node new?
          return;
        end
      end
      this.ports{end+1} = port;
    end

    function e = add_element(this,e)
      assert(~this.finalized); 
      for i=1:length(this.elements) % Is the element new?
        if(this.elements{i}==e)
          return;
        end
      end
      this.elements{end+1} = e;
    end

    function this = connect(this,port1,port2,varargin)
      assert(~this.finalized); 
      this.conns{end+1} = [{port1,port2},varargin];
    end

    % This function assign node ID, comports 'nodeNum' and 'maps'. 
    % The circuit node ID are assign as
    %   1. Assign ID circuit ports, by the order of add_port. 
    %      Note connected element ports have same ID with circuit ports. 
    %   2. Assign ID to visible ports of each element
    %   3. Assign ID to invisible nodes of each element. 
    function this = finalize(this)
      if(this.finalized), return; end

      % 1. Find all circuit ports and element ports
      nodes = this.ports;                % circuit ports
      for i=1:length(this.elements)      % element ports
        p = this.elements{i}.ports; 
        nodes(end+1:end+length(p)) = p;
      end
      M = length(nodes); inds = (1:M)';  % assign 1:M to all ports initially 
      %
      % 2. Find connected ports by 'conns'
      same = cell(length(this.conns),1); % nodes in each "conns" are the same 
      for i=1:length(this.conns)
        cports = this.conns{i};          
        cid = zeros(length(cports),1);
        for j=1:length(cports)
          for k = 1:M                    % find ID for nodes used in 'conn'
            if(cports{j}==nodes{k})      
              cid(j) = k; break;
            end
          end
        end
        if(any(cid==0)), error('ports not found'); end
        same{i} = cid;                   % all IDs should be merged
      end 
      %
      % 3. Merge IDs for connected ports  by 'same'
      allsame = false;                   % all 'conn' info applied?
      while(~allsame) 
        allsame = true;
        for i=1:length(same)
          cid = same{i}; 
          if(length(unique(inds(cid)))~=1) % not merged yet
            inds(cid) = min(inds(cid));  % all nodes use the min ID 
            allsame = false;
          end % if
        end % for
      end % while
      % 
      % 4. Reassign ID. Connected element ports use same ID with circuit ports, 
      %    followed by un-connected element ports
      [uinds,~,J] = unique(inds);        % inds = uinds(J) 
      N = length(uinds); inds = J;        % J is the ID assignment (M -> N) 
      %
      % 5. Compute 'nodeNum' and 'maps'. Internal nodes are assigned after N with element order 
      S = length(this.ports);            % start position of element ports in 'inds'
      maps = cell(length(this.elements),1);
      for i=1:length(this.elements)
        e = this.elements{i}; 
        np = length(e.ports); ni = e.nodeNum - np; assert(ni>=0);
        pid = inds(S+(1:np));            % element ports assigned as above 
        iid = N+(1:ni)';                 % element internals are appended at the end 
        maps{i} = [pid;iid];             % column vector
        N = N+ni; S = S+np;              % add number of internal nodes 
      end

      this.nodeNum = N; this.maps = maps;
      this.finalized = true; this.conns = {}; % clear temp data
      %if(this.is_leaf), this.flattened = true; end % opt for leaf?
    end % function

  end % method

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % I. Public methods for simulation interface 
  %   NOTE: all public methods support varargin to be used by subclass. 
  %
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Query interface. 
    %   For leaf circuits, return false by default. Must be overridden if provided .
    function is = ifc_simu(this) % support simulation?
      assert(this.finalized);
      if(this.is_leaf), is = false;   
      else
        is = true;
        for i=1:length(this.elements)
          if(~this.elements{i}.ifc_simu)
            is = false; return;
          end
        end
      end
    end

    function is = is_vsrc(this) % is it voltage source
      assert(this.finalized);
      if(this.is_leaf), is = false;
      else
        is = true;
        for i=1:length(this.elements)
          if(~this.elements{i}.is_vsrc)
            is = false; return;
          end
        end
      end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Current function: 
    %   i = I(v,varargin)
    %     v:  node voltages. nx1 vecotor. Could be nxm matrix, where each col for a point 
    %     i:  node currents, same size with v 
    %
    % Implementation of I: 
    %   I  ==  I_seq by default 
    %      ==  I_vec if vectorized 
    %   I_seq:  call each element's I function.   
    %   I_vec:  group elements by class, call one function for all elements with same class.
    %             The circuit must be vectorized to use this. It's for better performance
    %   No difference for leaf-circuits 
    %
    function i = I(this,v,varargin)
      assert(~this.is_leaf); % leaf circuit must override I()
      if(~this.vectorized)
        i = this.I_seq(v,varargin{:});
      else
        i = this.I_vec(v,varargin{:});
      end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Capacitance function: 
    %   c = C(v,varargin)
    % Similar with I 
    %
    function c = C(this,v,varargin)
      if(nargin<2||isempty(v))
        v = zeros(this.nodeNum,1);
      end
      assert(~this.is_leaf); % leaf circuit must override C()
      if(~this.vectorized)
        c = this.C_seq(v,varargin{:});
      else
        c = this.C_vec(v,varargin{:});
      end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ODE models: dV 
    %   dv = vdot(v,varargin) 
    %   
    function dv = dV(this,v,varargin)
      assert(this.finalized);
      i = this.I(v,varargin{:}); c = this.C(v,varargin{:});
      % NOTE: cap could be zero for input nodes
      %assert(~any(cap(:)==0)); % ODE -> DAE
      dv = i./c;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Voltage source: V(t) 
    %   v = V(t)
    function v = V(this,t,varargin)
      assert(~this.is_leaf);  % leaf circuit must override V(t)
      assert(this.finalized);
      error('input sources must implement this function');
    end
  end % method for simulation interface

  methods(Access=private,Sealed=true)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Private functions for I/C 
    %
    % Compute currents by recursively call elements' I function
    function i = I_seq(this,v,varargin) 
      assert(this.finalized);
      if(this.is_leaf)  % same with I() for leaf circuits
        i = this.I(v,varargin{:}); 
        return;
      end
      i = zeros(size(v));
      for ind=1:length(this.elements)
        e = this.elements{ind}; assert(~e.is_vsrc); % TODO: not support voltage source
        m = this.maps{ind};
        ii = e.I(v(m,:),varargin{:});
        i(m,:) = i(m,:)+ii; % assume no voltage source
      end
    end

    % Merge function calls for the same class, replaced by Static function call
    function i = I_vec(this,v,varargin)
      if(this.is_leaf)                        % same with I() for leaf circuits
        i = this.I(v,varargin{:}); return
      end
      if(~this.vectorized),this.vectorize;end % make sure it's vectorized
      [N,P] = size(v); assert(N==this.nodeNum); % N: nodes, P: points
      i = zeros(size(v)); 
      for cid=1:max(this.grp_cids)
        E = this.elements(this.grp_cids==cid); % K element with same class
        M = this.grp_maps{cid}; [Ne,K] = size(M); % M is Ne x K
        V = zeros(Ne,P,K);                    % Ne x P x K
        for ind=1:K
          V(:,:,ind) = v(M(:,ind),:);
        end
        Ie = E{1}.I_objs(E,V,varargin{:});    % call static method by any object
        Ifull = zeros(N,P,K);                 % Ie is Ne x P x K, Ifull is N x P x K
        for ind=1:K
          Ifull(M(:,ind),:,ind) = Ie(:,:,ind); 
        end
        i = i+sum(Ifull,3);
      end
    end

    function c = C_seq(this,v,varargin)
      assert(this.finalized);
      if(nargin<2||isempty(v))
        v = zeros(this.nodeNum,1);
      end
      if(this.is_leaf) % same with C() for leaf circuits
        c = this.C(v,varargin{:});
        return; 
      end
      c = zeros(size(v));
      for i=1:length(this.elements)
        e = this.elements{i}; assert(~e.is_vsrc); % TODO: not support voltage source
        m = this.maps{i};
        cc = e.C(v(m,:),varargin{:});
        c(m,:) = c(m,:)+cc; 
      end
    end

    function c = C_vec(this,v,varargin)
      if(this.is_leaf)                        % same with C() for leaf circuits.
        c = this.C(v,varargin{:}); return
      end
      if(~this.vectorized),this.vectorize;end % make sure it's vectorized
      [N,P] = size(v); assert(N==this.nodeNum);
      c = zeros(size(v)); 
      for cid=1:max(this.grp_cids)
        E = this.elements(this.grp_cids==cid); % K element with same class
        M = this.grp_maps{cid}; [Ne,K] = size(M); % M is Ne x K
        V = zeros(Ne,P,K);                    % Ne x P x K
        for i=1:K
          V(:,:,i) = v(M(:,i),:);
        end
        Ce = E{1}.C_objs(E,V,varargin{:});    % call static method by any object
        Cfull = zeros(N,P,K);                 % I is Ne x P x K, Ifull is N x P x K
        for i=1:K
          Cfull(M(:,i),:,i) = Ce(:,:,i); 
        end
        c = c+sum(Cfull,3);
      end
    end
  end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Static function for vectorized computation
    %   I = I_objs(Objs,V,varargin)
    %   C = C_objs(Objs,V,varargin)
    %     I/C/V is of size N*P*K N: # of nodes; P: # of points; K: # of objs 
    % NOTE: the default solution uses for loop, so essentially I_vec==I_seq. 
    %   Leaf circuit must override these functions for better performance. 
  methods(Static)
    function I = I_objs(Objs,V,varargin)
      assert(length(Objs)==size(V,3))
      I = ones(size(V));
      for i=1:length(Objs)
        assert(~Objs{i}.is_vsrc); % TODO not support voltage source
        I(:,:,i) = Objs{i}.I_vec(V(:,:,i),varargin{:});
      end
    end

    function C = C_objs(Objs,V,varargin)
      assert(length(Objs)==size(V,3))
      C = zeros(size(V));
      for i=1:length(Objs)
        assert(~Objs{i}.is_vsrc); % TODO not support voltage source
        C(:,:,i) = Objs{i}.C_vec(V(:,:,i),varargin{:});
      end
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % II. Public methods for verification interface 
  %   NOTE: We assume capacitance is independent of voltage here .
  %
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Query function 
    %
    function is = ifc_verify(this)
      assert(this.finalized);
      if(this.is_leaf), is = false;
      else
        is = true;
        for i=1:length(this.elements)
          if(~this.elements{i}.ifc_verify)
            is = false; return;
          end
        end
      end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Linear differential inclusion(LDI) model for currents
    %   [c,err] = I_ldi(region,vars,varargin)
    %     region:  specify the region for linearization, an object of Shape class 
    %     vars:    the index of nodes to be computed. 1:nodeNum by default 
    %     c,err:   LDI model: i = c*[v;1] +/- err
    % NOTE: We could use vars to speedup computation. 
    % TODO: remove vars from interface?
    %
    function [c,err] = I_ldi(this,region,vars,varargin)
      assert(~this.is_leaf); assert(this.finalized);
      n = this.nodeNum;
      if(nargin<3||isempty(vars)) % all nodes by default
        vars = 1:n; 
      end 
      vars = reshape(vars,1,[]); nv = length(vars);
      c = zeros(n,n+1); err = zeros(n,1);

      for i=1:length(this.elements)
        e = this.elements{i}; m = this.maps{i};
        ind = find(any(repmat(m,1,nv)==repmat(vars,length(m),1),2));
        if(~isempty(ind))
          rr = region.project(m); mm = m(ind);
          [cc,ee] = e.I_ldi(rr,ind,varargin{:});
          c(mm,[m;end]) = c(mm,[m;end])+cc;
          err(mm) = err(mm)+ee;
        end
      end

      c = c(vars,:); err = err(vars,:);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Linear differential inclusion(LDI) model for voltage
    %   [c,err] = dV_ldi(region,vars,varargin)
    %
    function [c,err] = dV_ldi(this,region,vars,varargin)
      assert(this.finalized);
      n = this.nodeNum;
      if(nargin<3||isempty(vars))
        vars = 1:n;
      end
      vars = reshape(vars,1,[]);
      [c,err] = this.I_ldi(region,vars,varargin{:});
      cap = this.C(zeros(n,1),varargin{:}); cap = cap(vars);
      c = c./repmat(cap,1,n+1); err = err./cap;
    end

  end % methods for verification interface


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % III.Public methods for small signal analysis 
  %
  methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Query function 
    %
    function is = ifc_ssan(this)
      assert(this.finalized);
      if(this.is_leaf), is = false;
      else
        is = true;
        for i=1:length(this.elements)
          if(~this.elements{i}.ifc_ssan)
            is = false; return;
          end
        end
      end
    end  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Jacobian matrix: 
    %   j = Jac(v,varargin) 
    %     v:  node voltage, nx1 vector. 
    %     j:  Jacobian matrix, nxn matrix.
    %
    function jac = Jac(this,v,varargin)
      assert(this.finalized);
      assert(size(v,2)==1);
      didv = this.dIdV(v,varargin{:}); % nxn
      c    = this.C(v,varargin{:}); % nx1
      jac  = didv/repmat(c,1,this.nodeNum);
    end

    % didv = dIdV(v,varargin) 
    function di = dIdV(this,v,varargin)
      assert(~this.is_leaf); assert(this.finalized);
      assert(size(v,2)==1);
      di = zeros(this.nodeNum);
      for i=1:length(this.elements)
        e = this.elements{i};
        m = this.maps{i};
        jj = e.dIdV(v(m,:),varargin{:});
        di(m,m) = di(m,m)+jj;
      end
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % IV. Functions that modify circuit tree 
  %   1. Vectorize   this = this.vectorize
  %   2. Flatten     this = this.flatten
  % NOTE: These functions change the circuit, please use copy to save the original one. 
  % TODO: modify or return a new circuit? It's more reasonable to modify 
  %
  methods(Access=public,Sealed=true)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Vectorize 
    %
    function this = vectorize(this) 
      assert(this.finalized)
      if(this.vectorized),return;end
      if(this.is_leaf)               % nothing to group
        this.vectorized=true; return
      end
      % NOTE: elements also modified, but elements are not shared (created in subclass constructor). 
      for i=1:length(this.elements)  % vectorize sub-circuits first
        this.elements{i}.vectorize;  % TODO: do or not? It's more reasonable to do
      end

      % 'grp_cids': assign class ID to each element 
      strs = {};                      % list of class name
      for i=1:length(this.elements)
        eid = this.elements{i}.vectorize_unique_id();
        inds = strcmp(eid,strs);      % find ID by class name
        if(any(inds))
          assert(length(find(inds))==1);
          this.grp_cids(i)=find(inds,1);
        else
          this.grp_cids(i)=length(inds)+1;
          strs{end+1} = eid;
        end % if
      end % for

      % 'grp_maps': calculate group mapping 
      for cid=1:max(this.grp_cids)
        eIDs = find(this.grp_cids==cid); assert(~isempty(eIDs));
        E = this.elements(eIDs);           % all elements of same class
        M = zeros(E{1}.nodeNum,length(E)); % col for each elements
        for i=1:length(eIDs) 
          M(:,i)=this.maps{eIDs(i)};       % col for each element
        end
        this.grp_maps{cid} = M;
      end

      this.vectorized=true;
    end % function

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Flatten 
    %
    function this = flatten(this) 
      assert(this.finalized);
      if(this.flattened), return; end;
      if(this.is_leaf)                    % nothing to do for leaf circuit 
        this.flattened=true; return; 
      end       
      if(this.vectorized)
        warning('Flatten destroys grouping information, please do vectorization again'); 
      end  
      while(~this.flattened)
        this.flatten_by_one_level; 
      end
      this.name = [this.name,'_flattened'];  % modify the circuit name
    end
  end
  
  %
  methods(Access=protected,Sealed=false)
    % subClass should override this function if it support different vectorize function with different configurations
    function id = vectorize_subtype_id(this)
      id = [];
    end
  end % methods

  % Private method for vectorize & flatten.
  %
  methods(Access=private,Sealed=true)
    % This function determine a unique ID of ojects for grouping during vectorization. 
    % The ID is defined to be className+subTypeName, where className must be unique for each class.
    % The subTypeName is use for distinguishing objects of the same class but different configurations.
    function id = vectorize_unique_id(this)
        id = class(this); 
        subid = this.vectorize_subtype_id;
        if(~isempty(subid))
          id = [id,'_',subid];
        end
    end
    %
    % This function flattens the circuit hierarchy by 1 level for the flatten method.
    % Comparing the new and original circuit, we have 
    %        orignal -> new
    %   name       same (or change?)
    %   nodeNum    same 
    %   ports      same (or not? same is better for identical interface)
    %   elements   elements of all elements
    %   maps       change*
    % maps{i} keeps the index mapping of the i-th elements and the orignal circuit, i.e.
    %   M1(i): Ei -> C
    % Ri.maps{j} keeps the mapping of the Ei's j-th elements and Ei, i.e. 
    %   M2(j): EEj -> Ei
    % The new maps{i,j} for new circuit should map the ij-th elements to the circuit nodes, i.e.
    %   M(i,j): EEj -> C == M1(M2)
    %
    function this = flatten_by_one_level(this)
      assert(this.finalized);
      % check if all elements are leaf circuits
      need_flatten = false;
      for i=1:length(this.elements)
        if(~this.elements{i}.is_leaf)
          need_flatten = true; break;
        end
      end

      % if all are leaf circuits, do not change.
      if(~need_flatten)
        this.flattened=true; return; 
      end
      
      % update elements
      new_eles = {};
      for i=1:length(this.elements)
        e = this.elements{i};
        if(e.is_leaf)
          new_eles{end+1} = e;
        else
          new_eles(end+1:end+length(e.elements)) = e.elements; 
        end
      end

      % update maps
      new_maps = {}; 
      for i=1:length(this.elements)
        Mi = this.maps{i}; e = this.elements{i};
        if(e.is_leaf)
          new_maps{end+1} = Mi;
        else % can be vectorized, but difficult to understand. who cares the flatten performance?
          for j=1:length(e.elements)
            Mj = e.maps{j}; new_maps{end+1} = Mi(Mj); 
          end
        end
      end

      this.elements = new_eles; this.maps = new_maps;
      this.vectorized = false; % need to redo grouping/vectorizing
    end
  end % methods


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % V. Utilities 
  %  TODO: get elements and ports?
  %    
  methods
    function isleaf = is_leaf(this)  % leaf circuit?
      assert(this.finalized);
      isleaf = isempty(this.elements); % without any subclass
    end

    function num = elemNum(this)      % number of elements
      num = length(this.elements);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % help for displaying internal structure
    %   index =   find_port_index(ports):      find port index
    %     ports:  an object of node class, could be a cell of nodes
    %     index:  the index of the port, 0 if not found
    %   names =   find_node_name(index);       return node names for given index
    %     index:  index of circuit nodes, all nodes by default 
    %     names:  a string if index is scaler, cell of strings if index is vector 
    %             node name if it's visible ports, 'internal_<index>' for internal nodes
    %   strs =    print_circuit_tree(prefix);  print the circuit tree structure   
    %     prefix: prefix of each string, '  ' by default
    %     strs:   a copy of strings printed.  
    %   strs =    print_status(prefix);        print circuit and sub-circuit status
    %     prefix: prefix of each string, '  ' by default
    %     strs:   a copy of strings printed.
    %   strs =    print_nodes;                 description of all nodes
    %     strs:   print the node name and their connections to sub-circuits 
    %
    function index = find_port_index(this,ports)
      assert(this.finalized);
      if(length(ports)==1&&~iscell(ports))  % convert to cell
        ports = {ports};
      end
      index = zeros(length(ports),1);
      for i=1:length(ports) 
        port = ports{i};
        for j=1:length(this.ports) 
          if(this.ports{j}==port) 
            index(i) = j; break; 
          end
        end
      end
    end

    function names = find_node_name(this,ind)
      if(nargin<2||isempty(ind)), ind=1:this.nodeNum; end
      if(numel(ind)==1&&ind<=0), ind = 1:this.nodeNum; end
      if(any(ind<=0)), error('index must be positive');end
      names = {};
      for i=1:length(ind)
        if(ind(i)<=length(this.ports))
          name = this.ports{ind(i)}.name;
        else
          name = ['internal_',num2str(ind(i)-length(this.ports))];
        end
        names{i} = name;
      end
      % return string for one and cell for a set
      if(length(names)==1), names = names{1}; end
    end

    function strs = print_circuit_tree(this,prefix)
      if(nargin<2),prefix='  '; end;
      % in-order tree traversal
      str = sprintf('%s%s@%s\n',prefix,this.name,class(this)); 
      strs{1} = str; fprintf(str); 
      for i=1:length(this.elements)
        estrs = this.elements{i}.print_circuit_tree([prefix,'  ']);
        strs(end+1:end+length(estrs)) = estrs;
      end
    end

    function strs = print_status(this,prefix)
      if(nargin<2),prefix='  '; end;
      str = sprintf('%s%s@%s:f=%d,v=%d\n',prefix,this.name,class(this),this.flattened,this.vectorized);
      strs{1} = str; fprintf(str);
      for i=1:length(this.elements)
        estrs = this.elements{i}.print_status([prefix,'  ']);
        strs(end+1:end+length(estrs)) = estrs;
      end
    end

    function strs = print_nodes(this)
      assert(this.finalized);
      % convert maps to matrix
      conn = zeros(this.nodeNum,length(this.elements));
      for i=1:length(this.maps)
        m = this.maps{i}; conn(m,i) = 1;
      end

      strs = {};
      str = sprintf('---- Circuit nodes of %s: ----\n',this.name);
      strs{1} = str; fprintf(str); 
      for i=1:this.nodeNum 
        % Visible ports
        if(i==1)
          str = sprintf('%d visible ports:\n',length(this.ports));
          strs{end+1} = str; fprintf(str);
        end
        str = sprintf('  %s\n',this.find_node_name(i)); 
        strs{end+1} = str; fprintf(str);
        es = find(conn(i,:)); % connected elements
        for j=1:length(es)
          eind = es(j); e = this.elements{eind}; m = this.maps{eind};
          epind = find(m==i);  % in case connected to multiple nodes
          for k=1:length(epind)
            str = sprintf('    <->%s of %s@%s\n',e.find_node_name(epind(k)),e.name,class(e)); 
            strs{end+1} = str; fprintf(str);
          end
        end
        if(i==length(this.ports)) 
          str = sprintf('%d internal nodes:\n',this.nodeNum-length(this.ports)); 
          strs{end+1} = str; fprintf(str);
        end
      end % for
    end % function

  end % methods

end % class
