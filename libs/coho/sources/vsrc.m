% Indepdent vtage source
%  vsrc(name,v)
%  v_p - v_m is always v
classdef vsrc < coho_vsrc 
  properties (GetAccess='private', SetAccess='private');
    vol;
  end
  methods
    function this = vsrc(name,v,ctg)
      if(nargin<2), error('not enough parameters'); end
      if(nargin<3||isempty(ctg)), ctg = true; end
      this = this@coho_vsrc(name,ctg);
      this.vol = v;
      this.finalize;
    end
    % v is v_p - v_m
    function v = V(this,t,varargin) 
      if(nargin<1||isempty(t)), t = 0; end;
      v = repmat(this.vol,1,numel(t)); 
    end
  end
end

