function err = lft_simu_err(model,bbox,c,lp)
% err = lft_simu_err(model,bbox,c,lp)
% This function computes bounds of linearization error
%   err(1) = min_{i IN bbox} model.data(i) - c'*v(i);
%   err(2) = max_{i IN bbox} model.data(i) - c'*v(i);
% Parameters
%   model: model loaded by ccm_getModel 
%   bbox:  the bounding box of regions
%    c:      coefficient
%   lp:     a linear program (optimal)
%   

if(nargin<4), lp = []; end;

% If bbox is small, lft_simu_err_bf is more efficient because the code is simple and vectorized.
% However, when bbox is large, lft_simu_err_bf is much slower because it uses large memory. 
% Mark developed a method which take advantages of convexity of nmos and pmos current function.
% Therefore, we switch to lft_simu_err_conv function for 'nmos' and 'pmos' models when bbox is large. 
% Generalize the function for other models if it is required. 

if( (isfield(model.META,'lib') && strcmpi(model.META.lib,'coho')) && ...                 % 'coho' library
    (isfield(model.META,'name') && any(strcmpi(model.META.name,{'nmos','pmos'}))) && ... % nmos or pmos models
    (~isfield(model,'err') || isempty(model.err) || all(model.err(:)==0)) && ...         % not interval based table
    (prod(max(1,diff(bbox,1,2)./model.GRID.dv)) > 2e5) )                                 % large region 
  err = lft_simu_err_conv(model,bbox,c); % ignore lp
else
  err = lft_simu_err_bf(model,bbox,c,lp);
end

