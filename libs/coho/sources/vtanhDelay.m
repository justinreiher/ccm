%  v = vtanhDelay(name,b,A,gain)
%    v = A*tanh(gain*(t-delay))+b
%    By default, delay is 0, b is 0, A is 1 and gain is 1e10; 
%    this voltage sources when evaluated takes as an input parameter the
%    delay, used in the nested bisection algorithm.
classdef vtanhDelay < coho_vsrc 
  properties (GetAccess='private', SetAccess='private');
    b,A,gain
  end
  methods
    function this = vtanhDelay(name,A,b,gain,ctg)
      if(nargin<2), error('not enough parameter'); end
      if(nargin<3||isempty(b)), b=0; end;
      if(nargin<4||isempty(gain)),gain = 5e10; end;
      if(nargin<5||isempty(ctg)), ctg=true; end;
      
       params = struct('amplitude',A,'offset',b,...
            'gain',gain);

      this = this@coho_vsrc(name,ctg,params);
      this.b = b; this.A = A;
      this.gain = gain;
      this.finalize;
    end
      
    % v_p - v_m
    function v = V(this,t,delay,varargin) % t can be a vector and so can delay, if t and delay are vectors they must be of the same dimension
      t = reshape(t,[],1);
      v = this.A*tanh(this.gain*(t-delay))+this.b; 
    end
    function params = getParams(this)
        params = struct('amplitude',this.A,'offset',this.b,...
            'gain',this.gain);
    end
  end
end
