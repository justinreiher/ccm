function [c,err] = lft_device(device,shape,wid,rlen,method,varargin)
% [c,err] = lft_device(device,shape,wid,rlen,method,varargin)
% This function computes a linear differential inclusion 
%   for a device function within a region specified by bbox and 
%   an optional linear program.
%   Y = c*X+[-err,err], s.t. X \in shape.
% Inputs:
%   device is the unique device name 
%   shape specify the region for computing LDI. It is a structure with two fields
%     'bbox':a nx2 matrix where bbox(:,1/2) is the lower/upper bound.
%     'lp':  an optional linear program of the form Ax <= b. 
%            It is recommended but not required to be coho LP.
%   wid: the width of the device, use default width if not provided
%   rlen: relative length of the device compared to the minimum lenght, '1' by default 
%        e.g. if process is 180nm and rlen is 2, then the device length is 180nm*2
%     NOTE: we assume the ids function is linear with the size of device.  
%   method is a string that specify the linfit algorithm.
%     for 'simu' models
%       'lls': linear least square method 
%       'mm':  find the minimum and maximum value (constant DI) 
%              error term is usually much larger.
%       if not provided, use ccm_cfg('get','lftSimuMethod');
%     for 'quad' models
%      'pt':   compute linear DI for quadratic model with optimal L2 norm error
%      'lip':  compute an approximated result based on Lipschitz constant
%              it is usually more efficient
%       if not provided, use ccm_cfg('get','lftQuadMethod');
%   varargin: more options for different linfit methods. 
%     see lft_simu, lft_quad for details
%
% Outputs:
%   c:   is a n+1x1 column vector
%   err: a non-negative value
if(nargin<2), error('not enough parameters'); end
if(nargin<3), wid = []; end % use default
if(nargin<4), rlen = []; end % use default
if(nargin<5), method=[]; end

model = ccm_getModel(device);
bbox = shape.bbox; lp = shape.lp;
switch(lower(model.META.type))
  case 'simu'
    if(isempty(method)), method=ccm_cfg('get','lftSimuMethod'); end
    [c,err] = lft_simu(model,bbox,method,lp,varargin{:});
  case 'quad'
    if(isempty(method)), method=ccm_cfg('get','lftQuadMethod'); end
    [c,err] = lft_quad(model,bbox,method,lp,varargin{:});
  otherwise
    error('not supported models');
end

% Fix linfit problem
%   The coefficient should be zero when the lower and upper bounds of
%   a variable are the same. However, the result from lft may not be the zero
%   because the bounding box is bloated outside in many methods or the result
%   may not find the optimal result (or one of optimal solutions).
%   This may cause problem when the variable is gnd or vdd. Therefore, we shift
%   the non-zero coefficient to the constant term.
eind = ( abs(bbox(:,1)-bbox(:,2))<=eps );
if(any(eind))
  xm = mean(bbox(eind,:),2);
  c(end) = c(end) + c(eind)'*xm;
  c(eind) = 0;
end
assert(~any(isnan(c)) && ~any(isnan(err)));

% Device size
if(isempty(wid)), wid = model.SIZE.wid; end
if(isempty(rlen)), rlen = 1; end
k = (wid/model.SIZE.wid)/rlen; 
c = c*k; err = err*k;

