% This class defines the 'NMOS' for COHO library
% The circuit has three nodes: s,g,d
% To create a NMOS, need to provide
%   1. name: circuit name
%   2. wid:  circuit width
%   2. ctg:  source connected to ground, false by default
%   4. rlen: (circuit length)/(minimum length), use 1 by default
% E.g. n = nmos('nmos',1e-5,0,true)
classdef nmos < coho_leaf
  properties (GetAccess='public', SetAccess='private');
    s; g; d;
  end
  methods
    function this = nmos(name,wid,ctg,rlen)
      if(nargin<2||isempty(name)||isempty(wid))
        error('must provide name and width'); 
      end
      if(nargin<3||isempty(ctg)), ctg = false; end
      if(nargin<4||isempty(rlen)), rlen = 1; end

      if(ctg) % s connected to ground
        subinfo.device='snmos';
        subinfo.I_factor=[0;-1]; 
        subinfo.C_factor=[1;1]; 
      else
        subinfo.device='nmos';
        subinfo.I_factor=[1;0;-1]; 
        subinfo.C_factor=[1;1;1]; 
      end
      this = this@coho_leaf(subinfo,name,wid,rlen); 
      if(~ctg)
        this.s = this.add_port(node('s'));
      end
      this.g = this.add_port(node('g'));
      this.d = this.add_port(node('d'));
      this.finalize;
    end
  end
end
