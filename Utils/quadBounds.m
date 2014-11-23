function [ymin, ymax] = quadBounds(As,x0s,rs,y0s)
% [ymin, ymax] = quadBounds(As,x0s,rs,y0s)
% The function compute the lower bound and upper bound of qudratic function [x;1]'*A*[x;1]
%   Find the x IN x0 +- r, that minimizes
%     u' * A * u
%   where u = [x; 1]
% 
%   A: n^2 x np matrix. each column for one polynomial
%   x0s: nxnp matrix. each column is a point. 
%   rs: nxnp matrix.  each column is the range for one point. e
%       each row is for one variable.
%   y0s: y value of x0s
%   ymin,ymax: 1xnp vector, lower/upper bound of qudratic function.
%
%   Note, if
%     y = sum_{i=0}^d sum_{j=i}^{d} b(i,j) * x_i * x_j
%   Then, b(i,j) is an element of an upper-triangular matrix.
%
%   Now, let
%     a(i,i) = b(i,i)/2
%     a(i,j) = b(i,j),  if j ~= i
%   We can now write
%     y = sum_{i=0}^d sum_{j=0}^{d} a(i,j) * x_i * x_j
%   where the a(i,j)'s are the elements of a symmetric matrix.
%   This is the matrix, A, that we expect as a parameter.
%   For example, if
%     y = 3*x_1^2 + x_1*x_2 - 5*x_1*x_3 + 2*x_2^2 - 6*x_2*x_3 - 7*x_3^2 +
%           9*x_1 - 3*x_2 + 8*x_3 - 12
%   Then we would have
%   A = [ 3.0,  0.5, -2.5,   4.5;
%         0.5,  2.0, -3.0,  -1.5;
%        -2.5, -3.0, -7.0,   4.0;
%         4.5, -1.5,  4.0, -12.0 ];
%  Note that grad y = 2*A*u
% 

% This function is optimized for linfitQuadFitLip function

[n,np] = size(x0s);
u0s = [ x0s; ones(1,np)]; % add the constant

if(nargin<4||isempty(y0s))
    dv = repmat(u0s,n+1,1).*reshape(repmat(reshape(u0s,[],1),1,n+1)',(n+1)^2,[]);
    y0s = sum(As.*dv,1);
end;

% compute derivative
dy0s = 2*reshape(sum(reshape(As.*repmat(u0s,n+1,1),n+1,[]),1),n+1,[]); % 2*A*u
dy0s = dy0s(1:n,:);

rAs = reshape(As,n+1,[]);
hs = reshape(sum(abs(rAs(1:n,:)),1),n+1,[]);
hs = 2*hs(1:n,:);
hs = hs.*rs; %maximum derivative change

dysmin = dy0s - hs;
dysmax = dy0s + hs;

% find the bound
dy = sum(rs.*max(abs(dysmin),abs(dysmax)),1); % maximum change of y
ymin = y0s - dy;
ymax = y0s + dy;

