% This class defines a capacitor with constant capacitance 
% It has two terminals 'x' and 'y'. 
% To create a capacitor, need to provide
%   1. name
%   2. c:  capacitance value
%   3. ctg: terminal 'y' connected to ground, true by default. 
% E.g. c = capacitor('cap',1e-12,true) 
%
classdef capacitor < circuit 
  properties (GetAccess='public', SetAccess='private');
    c; % capacitance
    x,y; % node 
  end
  methods
    function this = capacitor(name,c,ctg)
      if(nargin<3||isempty(ctg)), ctg=true; end
      this = this@circuit(name);
      this.c = c; assert(c~=0); % Support negative capacitor
      this.x = this.add_port(node('x'));
      if(~ctg) % connect to gnd otherwise
        this.y = this.add_port(node('y'));
      end
      this.finalize;
    end

    % support simulation
    function is = ifc_simu(this)
      is = true;
    end

    function c = C(this,v,varargin)
      c = this.c*ones(size(v));
    end

    % NOTE: for capacitor, I_c = c*dv/dt. 
    % However, our circuit class caculates dv/dt = I_other/c. 
    % By KCL, I_other + sum(I_c) = 0.
    % So new capacitor doesn't contribute to I_other.
    % Add a new capacitor changes dv/dt, but not I_other.
    function i = I(this,v,varargin)
      i = zeros(size(v));
    end

    % support small signal analysis 
    function is = ifc_ssan(this)
      is = true;
    end

    function didv = dIdV(this,v,varargin)
      didv = zeros(size(v));
    end

    % support verification
    function is = ifc_verify(this)
      is = true;
    end

    function [c,err] = I_ldi(this,region,varargin)
      c = zeros(this.nodeNum,2); err = zeros(this.nodeNum,1);
    end
  end

  methods(Access=protected)
    function id = vectorize_subtype_id(this)
      id = ['t_',num2str(this.nodeNum)]; 
    end
  end % methods

  methods(Static)
    function I = I_objs(Objs,V,varargin)
      I = zeros(size(V));
    end
    function C = C_objs(Objs,V,varargin)
      [N,P,K] = size(V); assert(length(Objs)==K);
      for i=1:K 
        C(1) = Objs{i}.c; 
      end
      C = repmat(C,[1,P,K]);
    end
  end
end
