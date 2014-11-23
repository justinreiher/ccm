function [c,err] = lft_simu_ls(model,bbox,lp)
% [c,err] = lft_simu_ls(model,bbox,lp)
%   This function compute the best linear fit using local search method.
%   Find the best linear fit of the ids function.
%       min err
%       err = max(abs(ids-c*[v;1]))
%   method:
%    1:  Efficient when region is small
%    2:  Efficient when region is large. However, lp is ignored and may not find the best fitting.
% NOTE:
% This function is expensive and the result is not good because of fmincon error
% The code does not work now because of fmincon error in the new Matlab version. 
% see http://www.mathworks.com/matlabcentral/newsreader/view_thread/174200 
%% Algorithm
%   For each grid point, error is convex with c, therefore, err(max/min
%   over all points) is convex with c. We can use fmincon to solve the problem.
%   We can always shift the constant term c(end) to make the err(1)=-err(2)
%   Therefore, the fmincon only search over a n-d space.

  opt = optimset('Display', 'off', 'LargeScale', 'off');
  scaleIds = 1e5; MIND=-0.1; MAXD=0.1; 
  
  n = length(model.GRID.v0);
  % find an initial point
  x0 = lft_simu_lls(model,bbox,lp);x0=x0(1:n);
  % opitmize L1 norm error
  f = @(c)(scaleIds*lft_simu_err_sym(model,bbox,[c;1],lp));
  c = fmincon(f,x0,[],[],[],[],MIND*ones(n,1),MAXD*ones(n,1),[],opt);
  % compute error 
  [err,c] = lft_simu_err_sym(model,bbox,[c;1],lp);

function [err,c] = lft_simu_err_sym(model,bbox,c,lp)
  if(nargin<4),lp=[];end
  err = lft_simu_err(model,bbox,c,lp);
  c(end) = c(end)+mean(err);
  err = diff(err)/2;
