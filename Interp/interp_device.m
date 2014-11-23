function Y = interp_device(device,X,wid,rlen,method)
% Y = interp_device(device,X,wid,rlen,method)
%   This computes the function value Y for a device of points X. 
% Inputs:
%   device: the unique name of the device 
%   X:      the points to be evaluated. dxn matrix, each column for a point
%   wid:    the width of the device, use the model size by default. 
%           can specify different width by nx1 vectors for each iterm of X. 
%   rlen:   relative length of the device compared to the minimum lenght, '1' by default 
%           e.g. for 180nm process, if rlen=2, then the device length is 2*180nm. 
%           can specify different length by nx1 vectors for each iterm of X. 
%   method: a string that specify the algorithm, see interp_data for details.
%           if method is not provided, ccm_cfg('get','interpIdsMethod') is used.
%   NOTE:   we assume the function is linear with the size of device.  
% Output:
%   Y:      function value, nx1 vector.
%
if(nargin<2), error('not enough parameters'); end
if(nargin<3), wid =[]; end
if(nargin<4), rlen =[]; end
if(nargin<5), method = ccm_cfg('get','interpIdsMethod'); end
[~,n] = size(X);
if(~any(numel(wid)==[0,1,n])), error('incorrect number of wid'); end
if(~any(numel(rlen)==[0,1,n])), error('incorrect number of rlen'); end

model = ccm_getModel(device);
switch(lower(model.META.type))
  case {'simu','quad'}
    Y = interp_data(model,X,method);
  otherwise
    error('not supported models');
end

if(isempty(wid)), wid = model.SIZE.wid; end
if(isempty(rlen)), rlen = 1; end
Y = Y.*(wid./model.SIZE.wid)./rlen; 
