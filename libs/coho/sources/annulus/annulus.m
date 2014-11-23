%   This class defines a brockett annulus for a signal. 
%       radius:     the radius of inner and outer ellipse.
%       vbnds:      The value of [V0l,V0h,V1l,V1h]. Usually, v0l =
%                   min(opoly(1,:)), v0h = min(ipoly(1,:)), v1l =
%                   max(ipoly(1,:)), v0h = min(opoly(1,:)). 
%       minT:       The minimal stay time of the stable state (1,3)
%       sbnds:      The brockett annulus is break into four discrete states
%                   1: low input 2: rising 3: high input 4: falling.
%                   sbnds(1) is the lower bound of state1, 
%                   sbnds(2) is the upper bound of state1, lower bound of 2 4
%                   sbnds(3) is the lower bound of state3, upper bound of 2 4
%                   sbnds(4) is the upper bound of state3
%                   It is common to set sbns = vbnds. 
%                   However, for COHO verification, the default value is
%                   sbnds(2) = vbnds(2)+0.05, sbnds(3) = vbnds(3)-0.05
%                   to make sure the linear approxiamted dv computed from 
%                   brockettfit function is clearly negative or positive.
%    polys:     The annulus is bounded by two ellipses. However, we use
%          polygons to approximate these ellipses for COHO verification.
%          It is optional. 
%                   
classdef annulus < coho_vsrc 
  properties (GetAccess='public', SetAccess='private');
    vbnds; sbnds;
    radius; polys;
    minT;
  end
  properties (GetAccess='private', SetAccess='private');
    trace;
  end

  methods
    % radius is required
    function this = annulus(name,radius,vbnds,minT,sbnds,polys) 
      if(nargin<2), error('not enough parameters'); end
      if(nargin<3) 
        vbnds =[0,0.15,1.65,1.8]; 
      end
      if(nargin<4||isempty(minT))
        minT = 0;
      end
      if(nargin<5||isempty(sbnds))
        sbnds = vbnds;
        sbnds(2) = sbnds(2)+0.05; % shift for linearization
        sbnds(3) = sbnds(3)-0.05;
      end 
      if(nargin<6||isempty(polys))  
        ll = vbnds(1); lh = vbnds(2); hl = vbnds(3); hh = vbnds(4); 
        load inner.mat; load outer.mat; 
        ipoly = [ lh + inner(1,:)*(hl-lh); inner(2,:)*radius(1) ]; 
        opoly = [ ll + outer(1,:)*(hh-ll); outer(2,:)*radius(2) ];
        polys = {ipoly,opoly};
      end

      this = this@coho_vsrc(name);

      this.radius = radius;
      this.vbnds = vbnds;
      this.minT = minT;
      this.sbnds = sbnds;
      this.polys = polys; % [ipoly,opoly]
      this.finalize;
    end

    % This function calculate the voltage at time t
    function v = V(this,t,varargin) 
      if(isempty(this.trace)) 
        error('call rand before simualtion'); 
      end 
      trace = this.trace; 
      ts = trace.ts; x0s = trace.x0s; as = trace.as; bs = trace.bs;

      % For the rise/fall phase, the function is 
      %   x(t) = (x0s(i)+bs(i)/as(i))*exp(as(i)*(t-ts(i)) - bs(i)/as(i) for t IN [ts(i),ts(i+1)]
      % For the stable phase, the function is 
      %   x(t) = x0s(i)+bs(i)*t for t IN [ts(i),ts(i+1)]
      n = length(t); v = zeros(n,1);
      for i=1:n % vectorize later
          ind = find(ts<=t(i),1,'last'); % find the interval
          if(isempty(ind)||ind==length(ts))
              error('simulation time is too long');
          end
          dt = t(i)-ts(ind); % relative time
          a = as(ind); b = bs(ind); x0 = x0s(ind);
          if(a~=0)
              v(i) = (x0+b/a)*exp(a*dt) - b/a; % rise/fall;
          else
              v(i) = x0+b*dt; % stable region
          end
      end
    end

    % This function compute LDI model for a Brockett annulus.
    % vbnds: the lower and upper bounds of current signal
    % region: the number of regions of current signal
    % method: linfit method to use
    function [c,err] = dV_ldi(this,vbnds,region,method)
      if(nargin<4), method = []; end; % use the default one
      [c,err] = lft_brock(this,vbnds,region,method);
    end

    % generate a random trace that satisfy the brockett annulus  
    %   T: the total time of the trace
    %   n: number of sampled points
    %  delay: delay of initial transistion  
    function this = rand(this,T,n,delay)
      if(nargin<3 || isempty(n))
          n = 20;
      end
      if(nargin<4 || isempty(delay))
          delay = 0;
      end
      
      ir = this.radius(1); or = this.radius(2);
      ll = this.vbnds(1); lh = this.vbnds(2); 
      hl = this.vbnds(3); hh = this.vbnds(4);
      minT = this.minT;
      assert((hh-ll)/n < (lh-ll));
      
      % the minimum transition time from lo/hi to hi/lo is 
      minTransT = (hh-ll)/(2*or)*(pi-2*acos((hl-lh)/(hh-ll)));
      hpt = (minT+minTransT);
      np = 2*ceil(T/(2*hpt));
      
      bnds = linspace(ll,hh,n+1)';
      xs = repmat(bnds(1:n),1,np) + repmat(diff(bnds),1,np).*rand(n,np);
      oys = sqrt( 1 - ((xs-(ll+hh)/2)/((hh-ll)/2)).^2 ) * or;
      iys = sqrt(max( 0, 1 - ((xs-(lh+hl)/2)/((hl-lh)/2)).^2 )) * ir; 
      ys = iys + (oys-iys).*rand(n,np);
      xs(:,2:2:end) =  xs(end:-1:1,2:2:end); % fall edge;
      ys(:,2:2:end) = -ys(end:-1:1,2:2:end);
      
      % transition
      x0 = xs(1:end-1,:); x1 = xs(2:end,:);
      y0 = ys(1:end-1,:); y1 = ys(2:end,:);
      a = (y1-y0)./(x1-x0);
      b = y0-a.*x0;
      t = (log(x1+b./a)-log(x0+b./a))./a;
      
      xend = ll+rand*(lh-hh);
      % stable
      st = (1+rand(1,np))*minT;
      slope = ([xs(1,2:end),xend] - xs(end,:))./st;
      as = [a;zeros(1,np)]; bs = [b;slope];
      x0s = xs; ts = [t;st];
      
      % initial 
      tinit  = delay*minT; xinit = ll+rand*(lh-ll); slope = (x0s(1)-xinit)/tinit;
      trace.as    = [0;reshape(as,[],1)];
      trace.bs    = [slope;reshape(bs,[],1)];
      trace.x0s   = [xinit;reshape(x0s,[],1);xend];
      trace.ts    = cumsum([0;tinit; reshape(ts,[],1)]);
      trace.phaseT = trace.ts(sort([2:n:end,1:n:end])); % [low rise high fall] order
      this.trace = trace;
    end

    % This function plots the annulus. 
    % showp: plot the polygons. 
    % showe: plot the ellipses.
    function display(this,showp,showe)
      if(nargin<2||isempty(showp))
        showp = true;
      end
      if(nargin<3||isempty(showe))
        showe = true;
      end

      
      hold on;
      if(showp) % plot the polygon
        ip = this.polys{1}; op = this.polys{2};
        plot(ip(1,[1:end,1]), ip(2,[1:end,1]),'k'); 
        plot(op(1,[1:end,1]), op(2,[1:end,1]),'k'); 
      end
      % plot the ellipse
      if(showe)
        % compute ellipse
        x = this.vbnds;
        ai = (x(3)-x(2))/2; 
        ao = (x(4)-x(1))/2;
        x0i = mean(x([2,3]));
        bi = this.radius(1);
        bo = this.radius(2); 
        x0o = mean(x([1,4]));
        xxi = x(2):0.001:x(3); 
        xxo = x(1):0.001:x(4); 
        yyi = annulus.ellipse_eval(ai,bi,x0i,xxi);
        yyo = annulus.ellipse_eval(ao,bo,x0o,xxo);
        plot(xxi,yyi,'b'); plot(xxi,-yyi,'b');
        plot(xxo,yyo,'b'); plot(xxo,-yyo,'b');
      end
    end

    % TODO do we need so much of parameters? 
    % clean up it latter. Check linfit/annulus
    function out = get(this,what)
      %   This function provide interface to access fields of brockett. Similar
      %   with get* function in OO programming. 
      %   'what' can be 
      %       'radius':   The raduis of inner/outer ellipse.
      %       'ipoly':    The inner bound polygon of the brockett
      %       'opoly':    The outer bound polygon of the brockett
      %       'iu/il/ou/ol''pl': The upper/lower polyline of the inner/outer polygon.
      %       'vbnds':    The voltage bounds of the inputs
      %       'v0/1l/h':  The value of lower/upper bound of 0/1
      %       'sbnds':    The stage bounds of the four discrete states
      %       'bnd12','bnd23','bnd34','bnd41':    The voltage bound of between states (1/2/3/4)
      %       'minT':     The minimal stay time of the stable states. 
      
      switch lower(what)
        % information of ellipse representation
        case 'radius'
          out = this.radius;
        case {'iradius','bi'}
          out = this.radius(1);
        case {'oradius','bo'}
          out = this.radius(2);      
        case 'x0i'
          out = mean(this.vbnds([2,3]));
        case 'x0o'
          out = mean(this.vbnds([1,4]));
        case 'ai'
          out = (this.vbnds(3)-this.vbnds(2))/2;
        case 'ao'
          out = (this.vbnds(4)-this.vbnds(1))/2;
      
        % information of polygon representation
        case 'ipoly'
          out = this.polys{1};
        case 'opoly'
          out = this.polys{2};
        case 'iupl' % split into polylines
          [~,out] = poly_split(this.polys{1},3);
        case 'ilpl'
          out = poly_split(this.polys{1},3);
        case 'oupl'
          [~,out] = poly_split(this.polys{2},3);
        case 'olpl'
          out = poly_split(this.polys{2},3);
      
        % information of regions
        case 'vbnds'
          out = this.vbnds;
        case 'v0l'
          out = this.vbnds(1);
        case 'v0h'
          out = this.vbnds(2);
        case 'v1l'
          out = this.vbnds(3);
        case 'v1h'
          out = this.vbnds(4);
        case 'sbnds'
          out = this.sbnds;
        case 'bnd12'
          out = this.sbnds(2);
        case 'bnd23'
          out = this.sbnds(3);
        case 'bnd34'
          out = this.sbnds(3);
        case 'bnd41'
          out = this.sbnds(2);
        case 'mint' % minimum stay time
          out = this.minT;
      
        % simulation trace
        case 'trace'
          out = this.trace;
          otherwise
              error([ what ' is not a parameter you can extract from a brockett' ]);
      end
    end
  end

  methods(Static) 
    function y = ellipse_eval(a,b,x0,x)
      y = (b/a).*(sqrt(a^2-(x-x0).^2));
    end
  end
end
