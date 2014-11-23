function [c, err]  = lft_simu_lls(model,bbox,lp)
% function [c, err] = lft_simu_lls(model,bbox)
%  Compute the linear approximation of the current which minimize the l2 error 
%  using linear least square method.
%  ids \in c'[v;1] +/- err, 
%  Parameters:
%    model:  model from ccm_getModel 
%       bbox specifies the domain cube for which we are working on:
%           bbox(i,1) is the lower bound for domain variable i.
%           bbox(i,2) is the upper bound for domain variable i.
%  Return:
%    c,err:  Linear approximation of the current of transistor
if(nargin<3),lp=[];end

% Find 'best' linear fit by least square method.
c = lft_simu_fit(model,bbox);

% Compute error term: err = actual - predicated
err = lft_simu_err(model,bbox,c,lp);

% shift the constant term
c(end) = c(end)+mean(err); err = diff(err)/2;
