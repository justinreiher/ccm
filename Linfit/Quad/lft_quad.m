function [c,err] = lft_quad(model,bbox,method,lp,useErr,varargin)
if(nargin<2)
  error('not enough parameter');
end
if(nargin<3 ||isempty(method)), method='lip'; end
if(nargin<4), lp = []; end
if(nargin<5), useErr = []; end

switch(lower(method))
  case 'pt' % brute force method
    [c,err] = lft_quad_pt(model,bbox,lp,useErr,varargin{:});
  case 'lip' % lipschitz approximation method
    [c,err] = lft_quad_lip(model,bbox,lp,useErr,varargin{:});
  otherwise
    error('not supported method\n');
end
