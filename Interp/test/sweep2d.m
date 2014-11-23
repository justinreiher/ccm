function [r,X1,X2] = sweep2d(model,V,method,t,grids,func)
% r = sweep2d(model,V,method,t,grids,func)
% This function is a wrapper function of interp_data/interp_jacob. 
% It sweep the voltage of two terminals specified by t and compute currents
% on these points.
%   model,V,method: the same with interp_data function. 
%   However, v is nx1 vector (one point)
%   t: indices of two terminals to be sweeped
%   grids: the sweep grids (a vector or  2x1 cell)
%   func: a function with interface func(model,V,method)
%   r: the return of the function
if(nargin<6||isempty(grids))
  eps=5e-2;
  grids = cell(2,1);
  for i=1:2 
    x0 = model.GRID.v0(t(i))+model.GRID.dv(t(i)); 
    x1 = model.GRID.v0(t(i))+model.GRID.nv(t(i))*model.GRID.dv(t(i)); 
    grids{i} = x0:eps:x1;
  end
end
if(~iscell(grids))
  grids = {grids,grids};  
end
if(nargin<7||isempty(func))
  func = @(model,V,method)interp_data(model,V,method);
end

n1 = length(grids{1}); n2 = length(grids{2}); n = n1*n2;
[X1,X2] = meshgrid(grids{1},grids{2});
V = repmat(V(:),1,n);
V(t,:) = [reshape(X1,1,[]);reshape(X2,1,[])]; 

r = func(model,V,method);
r = reshape(r,n1,n2);

% plot
if(nargout==0)
  figure; hold on; grid on;
  title([model.META.name,' Model']);
  xlabel(['v',num2str(t(1))]);ylabel(['v',num2str(t(2))]);zlabel('r');
  mesh(X1,X2,r); view(-37.5,30);
  axis([min(X1(:)),max(X1(:)),min(X2(:)),max(X2(:)),min(r(:)),max(r(:))]);
end
