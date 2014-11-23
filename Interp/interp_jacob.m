function [dYdX,Y] = interp_jacob(model,X)
% [dYdX,Y] = interp_jacob(model,X)
%   This function computes the Jacobian matrix of a function Y = f(X). i.e. dy/dx_i.
%   Here the function is implicit, expressed by sample points as in the model 
%
% Inputs:
%   model:   device model from ccm_loadModel.
%     model.GRID: the sampling grids.
%     model.data: the sampling data
%     model.META.type: model data type, 'simu' or 'quad'. 
%   X:       dxn matrix, each column for a value to be evaluated 
%
% Output:       
%    dYdX:  The Jacob matrix (dxn). Each column for a point 
%
% NOTE:
%   It uses 'coswin' method to compute ids(X) because others ('lookup','linear') are not C^1. 
% NOTE: 
%   It not recommend to use 'simu' models for computing Jacobian matrix. 
%   For the 'simu' models, the current function is a linear combination of sin/cos functions. 
%   (As current on grid points are constant, the jac is a1w'(x)+a2w'x(-1), where w is cos/sin functions).
%   For the 'simu' models, we can compute the Jacob matrix numerically using the "linear" interpolation. 
%   Although 'linear' interpolated current function is C^0, not C^1, it ususally works numerically.
%   See interp_jac_num

% check parameters 
[d,n] = size(X); 
if(d~=length(model.GRID.v0))
  error('incorrect dimension of model and voltage');
end

% compute the value of these point and grid information 
[Y,~,rems,gData,gY] = interp_data(model,X,'coswin');  

nc = 2^d; x = rems; xpi = x*pi; % dxnxnc 
% compute deriviatve of cosine window functions
cw = (128+150*cos(xpi)-25*cos(3*xpi)+3*cos(5*xpi))/256; % dxnxnc 
dcw= (-150*sin(xpi)+75*sin(3*xpi)-15*sin(5*xpi))*pi/256; % dxnxnc
pcw = zeros(d,n,nc);
for i=1:d
  ind = true(d,1); ind(i)=false;
  pcw(i,:,:) = prod(cw(ind,:,:),1);
end
% duplicate for d/dx_k (k=1:d)
w = repmat(prod(cw,1),[d,1,1]);
f = repmat(reshape(gY,[1,n,nc]),[d,1,1]);
% compute derivative of function f.
switch(model.META.type)
  case 'quad'
    % compute d(u'*A*u)/du = 2*u'*A for each point
    bs = reshape(cell2mat(gData),[(d+1)*(d+2)/2,n*nc]); % each column for a point
    xs = reshape(rems,[d,n*nc]); % position of each point
    df = reshape(evalQuadJac(bs,xs),[d,n,nc]);
  case 'simu'
    df = zeros(d,n,nc);
  otherwise
    error('Not supported type of models'); 
end
dYdX = pcw.*dcw.*f+w.*df;
dYdX= sum(dYdX,3); % dxn
dYdX= dYdX./repmat(model.GRID.dv,1,n);

