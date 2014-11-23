% ctg: m terminal connect to ground
classdef coho_vsrc < circuit
  properties (GetAccess='public', SetAccess='private');
    p,m;
  end
  methods
    function this = coho_vsrc(name,ctg)
      if(nargin<1||isempty(name)), error('not enough parameter'); end
      if(nargin<2||isempty(ctg)), ctg = true; end
      this = this@circuit(name);
      this.p = this.add_port(node('plus')); 
      if(~ctg)
        this.m = this.add_port(node('minus')); 
      end
    end
    function is = is_vsrc(this)
      is = true;
    end
    function is = ifc_simu(this)
      is = true;
    end

    % TODO: return 0 or NaN
    function i = I(this,v,varargin)
      i = zeros(size(v)); 
    end
    function c = C(this,v,varargin)
      if(nargin<2||isempty(v)), v = zeros(this.nodeNum,1); end
      c= zeros(size(v)); 
    end
  end % method
end % coho_vsrc

