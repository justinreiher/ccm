%  v = vsin(name,w,phi,b,A)
%    v = A*(sin(wt+phi)+b)
%    By default, phi is 0, b is 0, and A is 1;  
classdef vsin < coho_vsrc 
  properties (GetAccess='private', SetAccess='private');
    w,phi,b,A
  end
  methods
    function this = vsin(name,w,phi,b,A,ctg)
      if(nargin<2), error('not enough parameter'); end
      if(nargin<3||isempty(phi)), phi=0; end;
      if(nargin<4||isempty(b)), b=0; end;
      if(nargin<5||isempty(A)), A=1; end;
      if(nargin<6||isempty(ctg)), ctg=true; end;

      this = this@coho_vsrc(name,ctg);
      this.w = w; this.phi = phi; 
      this.b = b; this.A = A; 
      this.finalize;
    end
      
    % v_p - v_m
    function v = V(this,t,varargin) % t can be a vector
      t = reshape(t,[],1);
      v = this.A*(sin(this.w*t+this.phi))+this.b; 
    end
  end
end

