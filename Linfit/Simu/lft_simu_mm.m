function [c,err] = lft_simu_mm(model,bbox,lp)
% [c,err] = lft_simu_mm(model,bbox,lp,nodes) 
%  This function computes the max/min value of the ids (constant inclusion)
%  That is c(1:end-1) is zero, c(end) is the mean of max/min ids.
%    ids \IN c'*[x;1]+/-err  

n = length(model.GRID.v0);
c = zeros(n+1,1);
err = lft_simu_err(model,bbox,c,lp);
c(end) = c(end)+mean(err); err=diff(err)/2;
