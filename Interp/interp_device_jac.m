function dYdX = interp_device_jac(device,X,wid,rlen,method)
% dYdX = interp_device_jac(device,X,wid,rlen,method)
%   This computes dY/dX for the device function at points X
% Inputs:
%   device: the unique name of the device 
%   X:      the points to be evaluated. dxn matrix, each column for a point
%   wid:    the width of the device, use the model size by default. 
%           can specify different width by nx1 vectors for each iterm of X. 
%   rlen:   relative length of the device compared to the minimum lenght, '1' by default 
%           e.g. for 180nm process, if rlen=2, then the device length is 2*180nm. 
%           can specify different length by nx1 vectors for each iterm of X. 
%   method: a string that specify the algorithm
%           method == 'num', call interp_jacob_num function 
%           otherwise, call interp_jacob function 
%           if not provied, ccm_cfg('get','interpJacMethod') is used.
%   NOTE:   we assume the function is linear with the size of device.  
% Output:
%   Y:      function value, dxn vector.
%
if(nargin<2), error('not enough parameters'); end
if(nargin<3), wid =[]; end
if(nargin<4), rlen =[]; end
if(nargin<5), method = ccm_cfg('get','interpJacMethod'); end
[d,n] = size(X);
if(~any(numel(wid)==[0,1,n])), error('incorrect number of wid'); end
if(~any(numel(rlen)==[0,1,n])), error('incorrect number of len'); end

model = ccm_getModel(device);
switch(lower(model.META.type))
  case {'simu','quad'}
    if(strcmpi(method,'num')) 
      dYdX = interp_jacob_num(model,X); 
    else 
      dYdX = interp_jacob(model,X); 
    end
  otherwise
    error('not supported models');
end

% width and length 
if(isempty(wid)), wid = model.SIZE.wid; end
if(isempty(rlen)), rlen = 1; end
k = (wid./model.SIZE.wid)./rlen; 
if(length(k)==n), k = repmat(reshape(k,n,1),1,d)'; end
dYdX = dYdX.*k;  % dYdX is dxn matrix
