% This class defines a resitor with constant resistancee  
% It has two terminals 'x' and 'y'. 
% To create a capacitor, need to provide
%   1. name
%   2. r:  resistance value
%   3. ctg: terminal 'y' connected to ground, true by default. 
% E.g. r = resistor('resistor',1e-6,true) 
%
classdef resistor < circuit 
  properties (GetAccess='public', SetAccess='private');
    r; % resistance 
    x,y; % node
  end
  methods
    function this = resistor(name,r,ctg)
      if(nargin<2), error('not enough parameters'); end;
      if(nargin<3||isempty(ctg)), ctg = true; end;
      this = this@circuit(name);
      this.r = r; assert(r~=0); % Support negative resistor
      this.x = this.add_port(node('x'));
      if(~ctg)
        this.y = this.add_port(node('y'));
      end
      this.finalize;
    end

    % support simulation
    function is = ifc_simu(this)
      is = true;
    end

    % i_x = (v_y - v_x)/r
    function i = I(this,v,varargin)
      assert(size(v,1)==this.nodeNum); % two terminal device
      if(this.nodeNum==2) 
        i = [v(2,:)-v(1,:); v(1,:)-v(2,:)]./this.r;  % flow into is positive 
      else % connected to ground
        i = -v(1,:)./this.r;
      end
    end

    function c = C(this,v,varargin)
      c = zeros(size(v)); 
    end

    % support small signal analysis 
    function is = ifc_ssan(this)
      is = true;
    end

    % d(i_x)/d(v_x) = -1/r; d(i_x)/d(v_y) = 1/r; 
    function didv = dIdV(this,v,varargin)
      assert(size(v,1)==this.nodeNum); % two terminal device
      if(this.nodeNum==2)
        didv = [-1,1;1,-1]/this.r;
      else
        didv = -ones(size(v))./this.r;
      end
    end

    % support verification
    function is = ifc_verify(this)
      is = true;
    end
    
    function [c,err] = I_ldi(this,region,varargin)
      assert(size(v,1)==this.nodeNum); % two terminal device
      if(this.nodeNum==2)
        c = [-1,1; 1,-1]./this.r; err = zeros(2,1);
      else
        c = -1/this.r; err = zeros(1,1);
      end
    end
  end % methods

  methods(Access=protected)
    function id = vectorize_subtype_id(this)
      id = ['t_',num2str(this.nodeNum)]; 
    end
  end % methods

  methods(Static)
    function I = I_objs(Objs,V,varargin)
      [N,P,K] = size(V); assert(length(Objs)==K);
      rs = zeros(K,1);
      for i=1:K
        rs(i) = Objs{i}.r;
      end
      rs = repmat(reshape(rs,[1,1,K]),[N,P,1]);
      if(N==1)
        I = -V./rs; 
      elseif(N==2)
        I = [V(2,:,:)-V(1,:,:); V(1,:,:)-V(2,:,:)]./rs; 
      else
        error('incorrect input parameters');
      end
    end
    
    function C = C_objs(Objs,V,varargin)
      C = zeros(size(V));
    end
  end % methods
end
