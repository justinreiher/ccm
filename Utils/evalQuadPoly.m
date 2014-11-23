function v = evalQuadPoly(bs,xs)
% v = evalQuadPoly(bs,xs)
%  The function evaluate quadratic polynomial function 
%    v = b'*[x1^2;x1x2;...;x1xd;x1; x2^2;x2x3;...;x2xd;x2; ...; xd^2;xd; 1]
%      = [x;1]'*B*[x;1]
%    where B is an lower triangular matrix
%
% bs: each column is the coefficients for one polynomial function. (d+1)(d+2)/2 x n matrix
% xs: each column is a point. d x n matrix  
% v:  the value of each polynomial, nx1 vector.

[d,n] = size(xs);
if(size(bs,2)~=n || ~any(size(bs,1)~=[(d+1)*(d+2)/2,(d+1)^2]))
  error('The size of bs and xs does not match');
end  
% v = u'*B*u = sum(sum(B.*(u*u')));
%   = b'*(tril(u*u'));
u = [xs;ones(1,n)]; % d+1 x n
u1 = repmat(u,d+1,1); % (d+1)^2 x n. 
u2 = reshape(repmat(reshape(u,[],1),1,d+1)',(d+1)^2,n);
dv = u1.*u2; % each column is u*u' for a point

if(size(bs,1)==(d+1)*(d+2)/2)
  ind = tril(true(d+1,d+1)); % index of u*u';
  dv = dv(ind,:);  % (d+1)*(d+2)/2 x n
end
% v = b'*dv
v = sum(bs.*dv,1)'; 
