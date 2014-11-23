function [c,err] = lft_brock_ellipse(brockett,xbnd,region)
% [c,err] = lft_brock_ellipse(brockett,xbnd,region)
% The function computes the linear approximation of region bounded by 
%   two ellipses.
% y IN c*[x;1]+[-err,err]

% get inner/outer ellipses
ai = brockett.get('ai'); ao = brockett.get('ao');
bi = brockett.get('bi'); bo = brockett.get('bo');
x0i = brockett.get('x0i'); x0o = brockett.get('x0o');

switch(region)
  case {1,3}
    c = [0,0]; 
    err = max(ellipse_eval(ao,bo,x0o,xbnd));
  case {2,4}
    % NOTE when x is in a small interval, the computation error is large. 
    % Therefore, we set c as zero which is not optimal but error is tiny.
    xl = xbnd(1); xh = xbnd(2); 
    if(xh-xl >= 1e-4);
      P = [(xh^3-xl^3)/3, (xh^2-xl^2)/2; (xh^2-xl^2)/2, xh-xl];
      xxho = xh-x0o; xxlo = xl-x0o;
      iyo = (bo/(2*ao))*((xxho*sqrt(ao^2-xxho^2)+ao^2*asin(xxho/ao)) - ...
          (xxlo*sqrt(ao^2-xxlo^2)+ao^2*asin(xxlo/ao)));
      iyxo = x0o*iyo - (bo/(3*ao))*((ao^2-xxho^2)^(3/2)-(ao^2-xxlo^2)^(3/2));
      xxhi = xh-x0i; xxli = xl-x0i;
      iyi = (bi/(2*ai))*((xxhi*sqrt(ai^2-xxhi^2)+ai^2*asin(xxhi/ai)) - ...
          (xxli*sqrt(ai^2-xxli^2)+ai^2*asin(xxli/ai)));
      iyxi = x0i*iyi - (bi/(3*ai))*((ai^2-xxhi^2)^(3/2)-(ai^2-xxli^2)^(3/2));
      q = [iyxo+iyxi;iyo+iyi]./2;
      c = P\q; a = c(1); b = c(2);

      xo = x0o - (ao^2*a)/sqrt(ao^2*a^2 + bo^2);
      xi = x0i - (ai^2*a)/sqrt(ai^2*a^2 + bi^2);
      pto = [xl,xh,xo];
      do = a*pto+b-ellipse_eval(ao,bo,x0o,pto);
      pti = [xl,xh,xi];
      di = a*pti+b-ellipse_eval(ai,bi,x0i,pti);
      d  = [di,do];
      
      % shift constant term
      err = [min(d);max(d)];
      b = b-mean(err); err = diff(err)/2;
      c = [a,b];
    else % interval is small, use the middle point only
      vo = ellipse_eval(ao,bo,x0o,xbnd);
      vi = ellipse_eval(ai,bi,x0i,xbnd);
      y = [vo(:);vi(:)]; y = [min(y),max(y)];  
      c = [0,mean(y)]; err = diff(y)/2;
    end    
  otherwise 
    error('region must be 1-4');
end
if(region==4)
  c = -c;
end

function y = ellipse_eval(a,b,x0,x)
  y = (b/a).*(sqrt(a^2-(x-x0).^2));
