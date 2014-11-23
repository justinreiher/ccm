function jac = evalQuadJac(bs,xs)
% jac = evalQuadJac(bs,xs)
%  The function evaluates the derivative of a quadratic polynomial function as 
%   d(u'*A*u)/du = 2*u'*A for each point
% bs: each column is the coefficients for one polynomial function. (d+1)(d+2)/2 x n matrix
% xs: each column is a point. d x n matrix  
% jac:  the derivative of each polynomial, dxnn matrix.
[d,n] = size(xs);
if(any(size(bs)~=[(d+1)*(d+2)/2,n]))
  error('The size of bs and xs does not match');
end  
As = quadConvert(bs,'b','A'); % b -> A
us = [xs;ones(1,n)];
us = repmat(us,d+1,1); 
jac = 2*As.*us; % (d+1)^2 x n
jac = reshape(jac,[d+1,d+1,n]);
jac = reshape(sum(jac,1),[d+1,n]);
jac = jac(1:d,:); % dxn
