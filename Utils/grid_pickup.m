function [I,V] = grid_pickup(GRID,bbox,lp)
% [I,V] = grid_pickup(GRID,bbox,lp)
% This function find all grids points which are inside the bounding box 
% and satisfiy the LP if provided.
% Simulation parameters
%   GRID:  grid information of the circuit model 
%   bbox:  bounding box of the region
%   lp:    A linear program of the form Ax <= b. (not required to be cohoLP) 
%   I:     The indices of these grid points (a column vector)
%   V:     The voltage of these grids points (a column for a point)
if(nargin<3), lp = []; end;
if(isempty(bbox)), error('the bounding box is empty'); end
tol = 1e-6;

% Find all points in the bounding box
n = length(GRID.v0); % number of dimesions of the table
ibnds = grid_v2loc(GRID,bbox,repmat((1:n)',1,2));
ibnds = [floor(ibnds(:,1)), ceil(ibnds(:,2))];
% round-off error may cause problems 1) ibnds are out of range
% 2) ibnds(:,1)~=ibnds(:,2) when bbox(:,1)==bbox(:,2)=grid points. 
if(any(ibnds(:,1)<1-tol) || any(ibnds(:,2)>GRID.nv+tol))
  error('Out of range');
end
ibnds = [max(1,ibnds(:,1)),min(GRID.nv,ibnds(:,2))];

% Get ids, err and voltage for all grid points
inds = cell(n,1); vs = cell(n,1);
for i=1:n
  inds{i} = ibnds(i,1):ibnds(i,2);
  vs{i} = grid_loc2v(GRID,inds{i},i);
end
np = prod(diff(ibnds,[],2)+1);
cV = cell(n,1); [cV{:}] = ndgrid(vs{:}); 
cI = cell(n,1); [cI{:}] = ndgrid(inds{:}); 
V = zeros(n,np); I = zeros(n,np);
for i=1:n
  I(i,:) = reshape(cI{i},1,[]);
  V(i,:) = reshape(cV{i},1,[]);
end

% trim further by LP
trim = ~isempty(lp);
if(trim&&max(np*n,np*length(lp.b))>1e7), trim = false; end; % out of memory
% If the number of points is huge, the computation of lp may cause 'out of memory' error. 
% We can solve the problem by working on low dimensional case, e.g., 
% working on 2d plane for 3d models. See below for details.

if(trim)
  lp = lp_bloat(lp,GRID.dv-tol); % bloat lp to consider all bbox that intersect with LP
  ind = all(lp.A*V <= repmat(lp.b,1,np),1);
  V = V(:,ind); I = I(:,ind);
end

% convert to 1d index
cI = mat2cell(I,ones(1,n),size(I,2));
I = sub2ind(GRID.nv',cI{:})';
