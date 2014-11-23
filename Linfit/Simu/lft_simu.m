function [c,err] = lft_simu(model,bbox,method,lp)
% [c,err] = lft_simu(model,bbox,method,lp)
% This function find a linear fitting of the current of device within the specified region
%   | ids(v) - c'*[v;1] | <= err,  for any v within the cube (bbox) and lp.  
% Inputs
%   model:  model from ccm_getModel function
%   bbox:  The bounding box for each node. bbox(:,1/2) is the lower/upper bound
%   method: the algorithm to use  
%    'lls':  linear least square method. Efficient with larger error. Default method
%    'mm':   find the max/min ids, c(1:3) is zero, c(4) is mean of max/min ids.
%    'ls*':  local search method using convexity. Smaller error. Best linear fitting achieved. Very expensive
%    'lp*':  linear programming method, slow. Best linear fitting achieved. Too expensive to work.
%    *: not recommended for production usage
%   lp:  A linear program that specifies the region of the form lp.A*v <= lp.b. not required to be COHO LP
%
% Output
%   c:  The linear coefficent found
%   err: The maximum absolute error

if(nargin<2)
  error('not enough parameter');
end
if(nargin<3||isempty(method))
  method = 'lls';
end
if(nargin<4||isempty(lp))
  lp = [];
end

switch(lower(method))
  case 'lls' % Search the fitting which gives minimum l_2 error by linear least square method. 
    [c,err] = lft_simu_lls(model,bbox,lp);
  case 'mm' % Find the max/min value which means the coefficient is zero.
    [c,err] = lft_simu_mm(model,bbox,lp);
  case 'ls' % Search the best fitting which gives minimum l_inf error by convexity
    warning('The computation is too expensive');
    [c,err] = lft_simu_ls(model,bbox,lp);
  case 'lp' % Find the best fitting which gives minimum l_inf error by solving a lp
    warning('The function is too expensive to work');
    [c,err] = lft_simu_lp(model,bbox,lp);
  otherwise
    error('not supported method\n');
end

