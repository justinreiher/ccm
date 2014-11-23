function v = grid_loc2v(GRID,loc,dim)
% v = grid_loc2v(GRID,loc,dim)
% This function translates GRID locations to real values
% Inputs:
%   GRID:   model GRID, structure with v0 and dv
%   loc:    location on the GRID, nxm matrix 
%   dim:    which dim to work on, scalar or nxm matrix, 1 by default 
% NOTE: dim doesn't matter if all dimensions have same GRIDs 

if(nargin<3||isempty(dim)), dim = 1; end
v0 = reshape(GRID.v0(dim),size(dim));
dv = reshape(GRID.dv(dim),size(dim));
v  = loc.*dv+v0;
