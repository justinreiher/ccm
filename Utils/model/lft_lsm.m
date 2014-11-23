function [b,err,rerr] = lft_lsm(v,X)
% [b,err,rerr] = lft_lsm(v,X)
% The function computes the least square fit of function
%   Y = f(x1,...,xn);
%   min sum(err^2) 
%     err = Y - Xb;
%
% It implements the same function as \ provides. 
% However, it requires the last column of x to be 1
% And the last term of b is shifted to make the error balanced.
%
% v is a Mx1 vector, x is a Mxn matrix.
% v: The value of function on n points, nx1 vector.
% X: The variable of function.  nxk matrix, each column for one dimension. 
% b: the best linear coefficient.
% err: the maximum error.
% rerr: the maximum relative error

% The function computes the best least square fit of n variables function
%  Y = f(x1,...xn);
%
% min sqrt(sum(err^2))
%   err = Y - Xb
% 
% The solution is
%   X'Xb = X'Y
% 
% In fact, matlab's \ function implements this algorithm;
% 
if(~all(X(:,end)==1))
  error(' the last column of X must be 1 for constant term');
end;

b = X\v;

% shift constant
diffs = v - X*b;
err = [min(diffs); max(diffs)];
b(end) = b(end)+mean(err); 
err = diff(err)/2;
if(nargout>2)
  diffs = v - X*b;
  rerr = max(abs(diffs)./abs(v));
end;
