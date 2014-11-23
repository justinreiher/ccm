function [c,err] = lft_simu_lp(model,bbox,lp)
% [c,err] = lft_simu_lp(model,bbox,lp)
%   This function use linear programming to find the best linear fit.
%     min err
%     s.t. (for each grid)
%       c'V - err <= ids;
%       -c'V - err <=-ids;
%   Use interior method rather than simplex for speedup. 
%   Ususally, this method is slow when there are many grid points. 
%
% The method usually does not work 
%  1) it is quite expensive (large number of grids)
%   2) error of lp solver is huge

[ids,V] = grid_pickup(model.GRID,bbox,lp);

% construct the lp.
[n,np] = size(V);
U = [V;ones(1,np)];
A = [U', -ones(np,1);...
    -U', -ones(np,1)  ];
b = [ids;-ids];
f = [zeros(n+1,1);1];

%  option=optimset('MaxIter',1e4);
x = linprog(f,A,b,[],[],[],[],[],[]); % not coho LP (try cplex?)
c = x(1:n+1); err = x(end);
