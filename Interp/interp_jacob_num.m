function dYdX = interp_jacob_num(model,X,method,delta) 
% dYdX = interp_jacob_num(model,X,method,delta) 
% This function computes the jacob matrix using numerical method
%   method: could be 'lookup','linear','coswin'. 
%
% NOTE: the result may be nonsense for other method becase it is not C^1 function
%   The result is usually bad for 'lookup'; is OK most time for 'linear'.
%   For 'coswin' method, it works well for 'quad' model. 
%   For 'simu' model, it shares the same problem with interp_jacob.
% 
% NOTE, ususally the numerical method works
%      lookup  linear coswin
% simu   x       v    x
% quad   x       v    v

if(nargin<3), method= ccm_cfg('get','interpJacNumMethod'); end 
if(nargin<4), delta = 1e-12; end

[d,n] = size(X); dYdX = zeros(d,n); 
Xminus = X-delta; Xplus = X+delta; 
x0 = repmat(model.GRID.v0+model.GRID.dv,1,n); 
x1 = repmat(model.GRID.v0+model.GRID.dv.*model.GRID.nv,1,n);
ind0 = Xminus < x0; ind1 = Xplus > x1; 
Xminus(ind0) = x0(ind0); Xplus(ind1) = x1(ind1); 
dX = ones(d,n)*2*delta; dX(ind0|ind1) = delta;
for i=1:d
  X1 = X; X1(i,:) = Xplus(i,:); Y1 = interp_data(model,X1,method);
  X2 = X; X2(i,:) = Xminus(i,:); Y2 = interp_data(model,X2,method);
  dYdX(i,:) = (Y1 - Y2)';
end;
dYdX = dYdX./dX;
