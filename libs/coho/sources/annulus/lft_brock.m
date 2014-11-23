function [c,err] = lft_brock(brockett,ibnds,region,method)
% [c,err] = lft_brock(brockett,ibnds,region,method)
% The function computes linear approximation of input signal specified by 
%   Brockett's Annulus. (see Component/@brockett)
% Inputs 
%   brockett:   Brockett annulus from @brockett, which is bounded by two 
%         ellipses. It has an over-approximated polygon representation
%   ibnds:     lower/upper bounds of input signal
%   region:   discrete region (region 1, 2, 3, 4)
%   method:   linearization method
%            'lls': least square methods
%           'lp': linear programming method (call model_beltfit)
%           'ellipse': least square method on ellipse representation
% Outputs
%   c:       1x2 row vector, linear coefficient and constant term
%   err:     error term
%   y IN c*[x;1]+[-err,err]
% 
% NOTE
%   'lp' method can find the minimized L1 error based on the polygon 
%   representation. On the other hand, 'lls' method is faster and optimizes
%   the L2 error. 'ellipse' method works on the more accurate ellipse 
%   representation.     
%   It is not recommended to use 'lp' method for region 1 or 3. 
%   'c' is not around zero because 'cplex' always return the 'bad' optimal
%   point when there are multiple ones.

if(nargin<4 || isempty(method))
  method = 'lls'; 
end

% check the region
if(ibnds(1)>=ibnds(2))
  error('ibnds(1) can not be greater than ibnds(2)');
end
x = brockett.get('sbnds');
switch(region)
  case 1
    bnd = x([1,2]);
  case 3
    bnd = x([3,4]);
  case {2,4}
    bnd = x([2,3]);
  otherwise 
    error('region must be 1-4');
end
% NOTE: we bloat the invariant a little to compute the intersection plane. 
% Therefore, it is possible ibnds is a little outside the region. 
tol = 1e-3;
invalid = ibnds(1) > bnd(2)+tol | ibnds(2) < bnd(1)-tol; 
if(invalid)
  error('Region and ibnds does not match. The length of ibnds is %d',diff(ibnds));
else
  ibnds = min(bnd(2),max(bnd(1),ibnds));  
end

switch(lower(method))
  case 'ellipse'
    [c,err] = lft_brock_ellipse(brockett,ibnds,region);
  case 'lls'
    [c,err] = lft_brock_lls(brockett,ibnds,region);
  case 'lp'
    [c,err] = lft_brock_lp(brockett,ibnds,region);
  otherwise
    error('unknown method');
end

assert(~any(isnan([c,err]))); assert(~any(isinf([c,err]))); 
