%  vsrc = vpulse(name,range,period,delay)
%    create a pulse source.
%      range:  must be a [lo,hi], where lo is the low voltage level and hi is the high voltage level. 
%      period: must be [rt,ht,ft,lt]. 
%      delay:  0 by default 
%
%    The output of the pulse source is lo for t < delay.  It then transitions
%    linearly to hi over the next rt time units.  It remains high for
%    ht time units, and then returns to lo over the next ft time units and remains 
%    as low for lt time. The repeat the period.

classdef vpulse < coho_vsrc 
  properties (GetAccess='private', SetAccess='private');
    pwl,period,delay;
  end
  methods
    function this = vpulse(name,r,p,d,ctg)  
      if(nargin<3), error('not enough parameter'); end;
      if(nargin<4||isempty(d)), d=0; end;
      if(nargin<5||isempty(ctg)), ctg=true; end;
      if(numel(r)~=2), error('r must be [low,high]'); end;
      if(numel(p)~=4), error('p must be [rt,ht,ft,lt]'); end;
      if(~all(p>0)), error('p must be positive'); end;
      
      lo = min(r); hi = max(r);
      rt = p(1); ht = p(2); ft = p(3); lt=p(4);
      % piece-wise linear: time,rate,value 
      pwl = [[0,         (hi-lo)/rt,  lo];  % rise 
             [rt,        0,           hi];  % high 
             [rt+ht,     (lo-hi)/ft,  hi];  % fall 
             [rt+ht+ft,  0,           lo]   % low 
         ];

      params = struct('pwl',pwl,'period',sum(p),'delay',d);

      this = this@coho_vsrc(name,ctg,params);

      this.pwl = pwl; this.period = sum(p); this.delay = d;
      this.finalize;
    end
      
    function v = V(this,t,varargin) 
      t = reshape(t,[],1); pwl = this.pwl;
      p = mod(t-this.delay,this.period); % relative position of the period
      s = sum(repmat(p,1,size(pwl,1)) >= repmat(pwl(:,1)',length(t),1), 2); % index of stage
      v = (t >= this.delay) .* ( (p-pwl(s,1)) .* pwl(s,2) + pwl(s,3) )+ ...
          (t < this.delay) .* pwl(1,end); % low
      v = reshape(v,[],1);
    end

  end
end

