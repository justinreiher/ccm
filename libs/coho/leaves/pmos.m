% This class defines the 'PMOS' for COHO library
% The circuit has three nodes: s,g,d
% To create a PMOS, need to provide
%   1. name: circuit name
%   2. wid:  circuit width
%   2. ctp:  source connected to vdd, false by default
%   3. rlen: (circuit length)/(minimum length), use 1 by default
% E.g. p = pmos('pmos',1e-5,0,true)
classdef pmos < coho_leaf
  properties (GetAccess='public', SetAccess='private');
    s; g; d;
  end
  methods
    function this = pmos(name,wid,ctp,rlen)
      if(nargin<2||isempty(name)||isempty(wid))
        error('must provide name and width'); 
      end
      if(nargin<3||isempty(ctp)), ctp = false; end
      if(nargin<4||isempty(rlen)), rlen = 1; end

      if(ctp) % s connected to vdd
        subinfo.device='spmos';
        subinfo.I_factor=[0;-1];
        subinfo.C_factor=[1;1];
      else
        subinfo.device='pmos';
        subinfo.I_factor=[1;0;-1];
        subinfo.C_factor=[1;1;1];
      end
      this = this@coho_leaf(subinfo,name,wid,rlen);
      if(~ctp)
        this.s = this.add_port(node('s'));
      end
      this.g = this.add_port(node('g'));
      this.d = this.add_port(node('d'));
      this.finalize;
    end
  end
end
