% This class impelments the 'shape' abstract class for COHO circuits
classdef coho_shape < shape
  properties (GetAccess='public', SetAccess='protected');
    bbox;
    lp;
  end
  methods
    function this = coho_shape(p,varargin)
      if(nargin==1) % coho_shape(projectagon)
        this.bbox = p.bbox;
        this.lp = p.lp;
      else % coho_shape(bbox,lp);
        this.bbox = p;
        this.lp = varargin{1};
      end
    end
    function p = project(this,index)
      bbox = this.bbox(index,:);
      if(length(index)==length(unique(index)))
        lp = lp_pickup(this.lp,index);
      else % don't know how to project to duplicated variables.
        lp = [];
      end
      p = coho_shape(bbox,lp);
    end
    function b = rec(this)
      b = this.bbox;
    end
    % change variable x(index) to '-x(index)' 
    function p = neg(this,index)
      if(nargin<2||isempty(index))
        index = 1:size(this.bbox,1);
      end
      bbox = this.bbox; lp = this.lp;
      % flip bbox
      bbox(index,:) = -bbox(index,[2,1]);
      % flip variable in lp.
      lp.A(:,index) = -lp.A(:,index);
      lp.Aeq(:,index) = -lp.Aeq(:,index);
      p = coho_shape(bbox,lp);
    end
  end
end

