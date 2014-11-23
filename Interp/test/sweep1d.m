function [r,grids] = sweep1d(model,V,method,t,grids,func)
% [r,grids] = sweep1d(model,V,method,t,grids)
% This function is a wrapper function of interp_data/interp_jacob. 
% It sweep the voltage of one terminal specified by t and compute currents
% on these points.
%   model,V,method: the same with interp_data function. 
%   However, v is nx1 vector (one point)
%   t: the index of terminals to be sweeped
%   grids: the sweep grids
%   func: a function with interface func(model,V,method)
%   r: the return of the function
if(nargin<6||isempty(grids))
  eps = 5e-4;
  x0 = model.GRID.v0(t)+model.GRID.dv(t); 
  x1 = model.GRID.v0(t)+model.GRID.nv(t)*model.GRID.dv(t);
  grids = x0:eps:x1;
end
if(nargin<7||isempty(func))
  func = @(model,V,method)interp_data(model,V,method);
end

V = repmat(V(:),1,length(grids)); V(t,:) = reshape(grids,1,[]);
r = func(model,V,method);

% plot
if(nargout==0)
  figure; hold on; grid on;
  title([model.META.name,' Model']); 
  xlabel(['v',num2str(t)]); ylabel('r');
  plot(grids,r);
end
