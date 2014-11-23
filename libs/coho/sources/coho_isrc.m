classdef coho_isrc < circuit
  properties (GetAccess='public', SetAccess='private');
    p,m;
  end
  methods
    function this = coho_isrc(name,ctg)
      if(nargin<1||isempty(name)), error('not enough parameter'); end
      if(nargin<2||isempty(ctg)), ctg = true; end
      this = this@circuit(name);
      this.p = this.add_port(node('plus')); 
      if(~ctg)
        this.m = this.add_port(node('minus')); 
      end
    end
    function is = ifc_simu(this)
      is = true;
    end
    function c = C(this,v,varargin)
      c = zeros(size(v));
    end
  end % method
end % coho_isrc

