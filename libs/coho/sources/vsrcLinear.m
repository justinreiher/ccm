% Indepdent vtage source
%  vsrc(name,v)
%  v_p - v_m is always v
classdef vsrcLinear < coho_vsrc 
  properties (GetAccess='private', SetAccess='private');
    slope;
  end
  methods
    function this = vsrcLinear(name,slope,ctg)
      if(nargin<2), error('not enough parameters'); end
      if(nargin<3||isempty(ctg)), ctg = true; end
      params = struct('slope',slope);
      this = this@coho_vsrc(name,ctg,params);
      this.slope = slope;
      this.finalize;
    end
    % v is v_p - v_m
    function v = V(this,t,varargin) 
      if(nargin<1||isempty(t)), t = 0; end;
      v = this.slope*t; 
    end
  end
end

