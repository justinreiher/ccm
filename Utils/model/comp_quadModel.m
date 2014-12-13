function Mq = comp_quadModel(M,GRID,path)
% This funvtion computes quadratic polynomial model based on the simulation model.
% The center of (i,j,...) cell is  (v0+i*dv,v0+j*dv,...), the length of the cell is 2*dv;
if(nargin<3||isempty(path)), path= '.'; end

% check parameters
if(~strcmpi(M.META.type,'simu'))
  error('The model must be of type simu');
end
% number of dimensions.
d = length(M.GRID.v0); 
% cell GRID information
v0 = GRID.v0; nv = GRID.nv; dv = GRID.dv;
if(numel(v0)==1), v0 = repmat(v0,d,1); end;
if(numel(nv)==1), nv = repmat(nv,d,1); end;
if(numel(dv)==1), dv = repmat(dv,d,1); end;
if( any(size(v0)~=[d,1]) || any(size(dv)~=[d,1]) || any(size(nv)~=[d,1]) )
  error('incorrect GRID information');
end
% number of grids in a cell
eps = 1e-12;
ng = dv./M.GRID.dv;
if(any(mod(ng,1)>eps & mod(ng,1)<1-eps))
  error('incorrect GRID information: nv/M.GRID.nv should be an integer');
else
  ng = round(ng);
end
% find the position of v0 in M. 
v0m = grid_v2loc(M.GRID,v0,(1:d)');
if(any(mod(v0m,1)>eps & mod(v0m,1)<1-eps))
  error('incorrect GRID information: v0 should be GRID of M');
else
  v0m = round(v0m);
end
GRID.v0 = v0; GRID.nv = nv; GRID.dv = dv;

% compute quadratic model for each hyper-rectangle
N = prod(nv); % number of cubes
data = cell(N,1); err = zeros(N,1); inds = cell(d,1);
% for each cube, normalize it to [-1,1]^d, and compute polynomial terms for GRID points.
X = quad_gridX(2*ng+1);
cube = cell(d,1); 
for i=1:N
  [inds{:}] = ind2sub(nv',i); % index of the cell
  lo = v0m+([inds{:}]'-1).*ng; hi = v0m+([inds{:}]'+1).*ng;
  for j=1:d
    cube{j} = lo(j):hi(j);
  end
  V = M.data(cube{:});
  [c,u] = lft_lsm(V(:),X);  % with equal error
  data{i} = c; err(i) = u;
  if(~isempty(M.err))
    E = M.err(cube{:});
    err(i) = err(i)+max(E(:));
  end
end

data = reshape(data,nv'); err = reshape(err,nv');

% save result
META = M.META; SIZE = M.SIZE;
META.type = 'quad';  % update type
file = [path,'/',META.name,'.mat'];  
save(file,'data','err','GRID','SIZE','META');
Mq = load(file);

function x = quad_gridX(siz)
% x = interp_genGrid(siz)
%  The function break a cube with [-1,1] as the lower/upper bound into
%  grids defined by siz. And compute quadratic polynomial for thess points, 
%   i.e. [x1^2,x1x2,...,x1xd,x1, x2^2,...x2xd,x2, ... ,xd^2,  ,1] 
%
siz = reshape(siz,1,[]);
n = length(siz); np = prod(siz);
V = ones(np,n+1);
for i=1:n
  SIZE = ones(1,n); SIZE(i) = siz(i);
  SIZE2 = siz; SIZE2(i) = 1;
  v  = linspace(-1,1,siz(i));
  v = reshape(repmat(reshape(v,SIZE),SIZE2),[],1);
  V(:,i) = v; % xi of all grids
end;

ind = 1;
x = zeros(np,(n+1)*(n+2)/2);
for i = 1:n+1
  for j=i:n+1
    x(:,ind) = V(:,i).*V(:,j);
    ind = ind+1;
  end;
end;
