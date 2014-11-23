function err = lft_simu_err_bf(M,bbox,c,lp)
if(nargin<4), lp = []; end

% get indices of grids points in this region
[I,V]= grid_pickup(M.GRID,bbox,lp); 
% compute error 
e = M.data(I) - [V;ones(1,size(V,2))]'*c; % (c'*[v;1])' 
% add M error for interval table 
if(~isempty(M.err)) 
  e = [e - M.err(I); e+M.err(I)]; 
end  
err = [min(e);max(e)];
