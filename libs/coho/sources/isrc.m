% Indepdent current source
%  isrc(name,i)
%  i_p - i_m is always i
classdef isrc < coho_isrc 
  properties (GetAccess='private', SetAccess='private');
    i;
  end
  methods
    function this = isrc(name,i,ctg)
      if(nargin<2), error('not enough parameters'); end;
      if(nargin<3||isempty(ctg)), ctg = true; end
      this = this@coho_isrc(name,ctg);
      this.i = i;
      this.finalize;
    end
    % NOTE: no such function interface
    % NOTE: Do not support isrc?
    function i = I(this,v,varargin) 
      i = repmat([1;-1],1,size(v,2)).*this.i;
      if(this.nodeNum==1), i = i(1,:); end
    end
  end
end

