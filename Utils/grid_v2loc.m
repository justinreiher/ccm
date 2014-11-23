function loc = grid_v2loc(GRID,v,dim)
% loc = grid_v2loc(GRID,v,dim)
% This function translates real value to the location on the GRIDs 
% Inputs:
%   GRID:   structure with v0 and dv
%   v:      voltage value, nxm matrix 
%   dim:    which dim to work on, scalar or nxm matrix, 1 by default 
% NOTE: dim doesn't matter if all dimensions have same GRIDs 

if(nargin<3||isempty(dim)), dim = 1; end

v0 = reshape(GRID.v0(dim),size(dim));
dv = reshape(GRID.dv(dim),size(dim));
loc = (v-v0)./dv;

