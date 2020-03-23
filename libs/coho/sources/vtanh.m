%  v = vtanhClock(name,w,phi,b,A,gain,delay)
%    v = A*tanh(gain*sin(wt+phi))+b
%    By default, phi is 0, b is 0, and A is 1;  
classdef vtanh < coho_vsrc 
  properties (GetAccess='private', SetAccess='private');
    b,A,gain,delay
  end
  methods
    function this = vtanh(name,A,b,gain,delay,ctg)
      if(nargin<2), error('not enough parameter'); end
      if(nargin<3||isempty(b)), b=0; end;
      if(nargin<4||isempty(gain)),gain = 3; end;
      if(nargin<5||isempty(delay)),delay=1e-16; end;
      if(nargin<6||isempty(ctg)), ctg=true; end;
      
      params = struct('amplitude',A,'offset',b,'delay',delay,...
            'gain',gain);

      this = this@coho_vsrc(name,ctg,params);
      this.b = b; this.A = A;
      this.gain = gain; this.delay = delay;
      this.finalize;
    end
      
    % v_p - v_m
    function v = V(this,t,varargin) % t can be a vector
      t = reshape(t,[],1);
      v = this.A*tanh(this.gain*(t-this.delay))+this.b; 
    end
    function params = getParams(this)
        params = struct('amplitude',this.A,'offset',this.b,'delay',this.delay,...
            'gain',this.gain);
    end
  end
end
