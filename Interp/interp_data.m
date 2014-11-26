function [Y,inds,rems,gData,gY,gW] = interp_data(model,X,method) 
% Y = interp_data(model,X,method) 
%   This evaluates a function y = f(x) by interpolation method. 
%   The implicit function 'f' is modeled by sampled data as specified by 'model'. 
%   It interpolats the sampled data and returns the function value f(x). 
%
% Inputs:
%   model:   device model from ccm_loadModel.
%     model.GRID: the sampling grids.
%     model.data: the sampling data
%     model.META.type: model data type, 'simu' or 'quad'. 
%   X:       dxn matrix, each column for a value to be evaluated 
%   method:  interpolation methods
%     'lookup':   find the nearest grid point and use its value
%     'linear':   find all neighbors and use bilinear interpolation.
%     'coswin':   find all neighbors and use 'cosine window' interpolation.
%
% Output:       
%  Y:    The function value of points X. It is reshaped to a nx1 vector 
%
% NOTE: The continuity of these interpolated functions are
%                   simu   quad
%          lookup   x      x
%          linear   C^0   C^0
%          coswin   C^n   C^n
%
% [Y,inds,rems,gData,gY,gW] = interp_data(model,X,method) 
%   For internal usage in inter_jacob.m
% Output
%   inds,rems: position of X on the grids
%   gData:     model data on grid points
%   gY:        function value of related grid points
%   gWs:       weight of these grid points

% check parameters 
if(nargin<3||isempty(method)), method = 'coswin'; end
[d,n] = size(X); 
if(d~=length(model.GRID.v0))
  error('incorrect dimension of model and voltage');
end

% compute position in the table
dim = repmat((1:d)',1,n);
grids = grid_v2loc(model.GRID,X,dim); % table index
% NOTE: Simulator may query value outside the boundary.
% If v is out of boundary, we use the value of the nearest cube. 
grids = max(1,min(grids,repmat(model.GRID.nv,1,n)));

% comupte position and related grids dxnxnc marix 
%   'inds': position of related grid
%   'rems': relative positive of the point to the grid
switch(lower(method))
  case {'linear','coswin'}
    % find all 'close' grid points and use the linear/cosine interpolation of their values. 
    % We say a grid point g is close to a point p if all(abs(g-p)<=1) 
    % 1. find lower and upper grids for each dimension 
    lInd = floor(grids); hInd = lInd+1;
    lRem = grids - lInd; hRem = grids - hInd; % lRem - hRem = 1
    % It is OK to use model.GRID.nv because the gW is 0 for both linear and coswin method 
    hInd = min(hInd,repmat(model.GRID.nv,1,n));
    Ind = [lInd,hInd]; Rem = [lRem,hRem];
    % 2. find all (2^d) 'close' grid points by combination
    cInd = mat2cell(Ind,ones(1,d),n*ones(1,2)); % dx2 cell, each cell for all points
    cRem = mat2cell(Rem,ones(1,d),n*ones(1,2)); 
    ind = utils_combineInd(d,2,true); 
    nc = 2^d; % number of combinations
    cInd = cInd(ind); cRem = cRem(ind); %d dx2^d cell, each column is for one combination
    % 3. convert to d x n x nc matrix
    inds = reshape(cell2mat(cInd),[d,n,nc]); 
    rems = reshape(cell2mat(cRem),[d,n,nc]); 
  case 'lookup'
    % Find the nearest grid point and use its value
    nc = 1; 
    inds = round(grids); rems = grids - inds; % dxn
    inds = reshape(inds,[d,n,nc]); %inds = permute(inds,[2,3,1]);
    rems = reshape(rems,[d,n,nc]); %rems = permute(rems,[2,3,1]);
  otherwise
    error('do not support method'); 
end

% get data from model
ii = reshape(mat2cell(inds,ones(1,d),n,nc),[d,1]); % 1xnxnc
ii = reshape(sub2ind(model.GRID.nv',ii{:}),[n,nc]); % nxnc
gData = model.data(ii); % nxnc

% evaluate Y
switch(lower(model.META.type))
  case 'simu' 
    gY = gData;
  case 'quad'
    bs = reshape(cell2mat(gData),[(d+1)*(d+2)/2,n*nc]); % each column for a point
    xs = reshape(rems,[d,n*nc]); % position of each point
    gY = reshape(evalQuadPoly(bs,xs),[n,nc]);
  otherwise 
    error('Not supported type of models'); 
end

% compute gW nxnc matrix
switch(lower(method))
  case 'coswin'
    xpi = rems*pi;
    cw = (128+150*cos(xpi)-25*cos(3*xpi)+3*cos(5*xpi))/256;
    gW = reshape(prod(cw,1),[n,nc]);
  case 'linear'
    gW = reshape(prod(1-abs(rems),1),[n,nc]);
  case 'lookup'
    gW = ones(n,nc);
  otherwise
    error('do not support method'); 
end
Y = sum(gY.*gW,2);
