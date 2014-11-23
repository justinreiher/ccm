function [c,err] = lft_brock_lp(brockett,ibnds,region)
% [c,err] = lft_brock_lp(brockett,ibnds,region)
% Find the optimal L1 norm error by solving LPs using 
% polygon representation
  switch(region)
    case {1,3}
      below = brockett.get('olpl');
      above = brockett.get('oupl');
    case 2 
      below = brockett.get('iupl');
      above = brockett.get('oupl'); 
    case 4
      below = brockett.get('olpl');
      above = brockett.get('ilpl'); 
    otherwise 
      error('region must be 1-4');
  end
  [c, err] = lft_beltfit(below, above, ibnds); 
end


function [c, err] = lft_beltfit(below, above, x)
% [c, err] = lft_beltfit(below, above, x)
% This function compute a linear approximation of 
%   a belt bounded by upper/lower polylines.
% Inputs
%   below: 2xn matrix with each column as a point. 
%     it is a sequence of points that define a piecewise-linear function
%       whose value is always less than or equal to the function value 
%       we are approximating.  
%   above: similar with below, the upper bound of the function. 
%   x: a row vector with two points.  The first point is the left
%       endpoint for the interval in which the approximation is to hold.
%       The second point is the right endpoint of this interval.
% Outputs: 
%   c:  The coefficients for the linear approximation.
%       f(x) IN c * x + [-err,err] 
%   err:  The maximum error.
%
% Restrictions: 
%   The points in below(:,1) must be in increasing order.  
%   This assumption is not checked in this implementation.  
%   Likewise for above(:,1).  
%   x(1) and x(2) must be in the range spanned by below(1,:)
%   and above(1,:).  An error is thrown if this does not hold.

  % Make sure that our range is contained in the span of 'below' and 'above'.
  if( (x(1) < below(1,1)) || (x(2) > below(1,end)) || ...  
    (x(1) < above(1,1)) || (x(2) > above(1,end)) )
      error('x is out of range');
  end
  
  % Find all points in the range of [x(1),x(2)]. 
  apts = polyline_crop(above,x); na = size(apts,2);
  bpts = polyline_crop(below,x); nb = size(bpts,2);

  % Now, it's a simple matter of linear programming to find the best fit.
  Aa = -[apts(1,:)', ones(na,1), ones(na,1)]; % axi+b+u > yi
  ba = -apts(2,:)';
  Ab = [bpts(1,:)', ones(nb,1), -ones(nb,1)]; % axi+b-u < yi
  bb = bpts(2,:)';
  A = [Aa;Ab]; b = [ba;bb];
  f = [0;0;1];
  
  % solve LP by matlab solver
  [c,err] = linprog(f,A,b); 
  c = c(1:end-1)';  % c(end)==err
end 
