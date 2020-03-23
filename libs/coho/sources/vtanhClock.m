%  v = vtanhClock(name,w,phi,b,A,gain,delay)
%    v = A*tanh(gain*sin(wt+phi))+b
%    By default, phi is 0, b is 0, and A is 1;  
classdef vtanhClock < coho_vsrc 
  properties (GetAccess='private', SetAccess='private');
    w,phi,b,A,gain,delay,holdHigh
  end
  methods
    function this = vtanhClock(name,w,phi,b,A,gain,delay,holdHigh,ctg)
      if(nargin<2), error('not enough parameter'); end
      if(nargin<3||isempty(phi)), phi=0; end;
      if(nargin<4||isempty(b)), b=0; end;
      if(nargin<5||isempty(A)), A=1; end;
      if(nargin<6||isempty(gain)),gain = 3; end;
      if(nargin<7||isempty(delay)),delay=0; end;
      if(nargin<8||isempty(ctg)), ctg=true; end;
      if(nargin<9||isempty(holdHigh)), holdHigh = false; end;
      
      params = struct('frequency',w,'phaseShift',phi,'offset',b,...
            'amplitude',A,'gain',gain,'delay',delay,'holdHighFlag',holdHigh);

      this = this@coho_vsrc(name,ctg,params);
      this.w = w; this.phi = phi; 
      this.b = b; this.A = A;
      this.gain = gain; this.delay = delay;
      this.holdHigh = holdHigh;
      this.finalize;
    end
      
    % v_p - v_m
    function v = V(this,t,varargin) % t can be a vector
      t = reshape(t,[],1);
      v = heaviside(t-this.delay).*(this.A*tanh(this.gain*sin(this.w*(t-this.delay)+this.phi))+this.b);
      if((t-this.delay)==0)
          t = t+1e-16;
      end
      if(this.holdHigh)
          v = v+heaviside(-t+this.delay);
      end
      
    end
    
    function params = getParams(this)
        params = struct('frequency',this.w,'phaseShift',this.phi,'offset',this.b,...
            'amplitude',this.A,'gain',this.gain,'delay',this.delay,'holdHighFlag',this.holdHigh);

    end
  end
end
