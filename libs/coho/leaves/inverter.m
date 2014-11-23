% This class defines the 'inverter' circuit for COHO library 
% The circuit has two nodes: i,o
% To create an inveter, need to provide
%   1. name: circuit name
%   2. wid:  circuit width 
%            wid(1) is the nmos width
%            wid(2) is the pmos width, use 2*wid(1) if not provided.
%   3. rlen: (circuit length)/(minimum length), use 1 by default
% E.g. inv = inveter('inv',[1e-5,2e-5],1)
%
% NOTE: the mat file required is ['inv_',num2str(wid(2)/wid(1))]
classdef inverter < coho_leaf
  properties (GetAccess='public', SetAccess='private');
    i; o;
  end
  methods
    % r:    the ratio of pmos/nmos. make sure inv_<r>.mat is available
    % size: the size of the nmos.
    function this = inverter(name,wid,rlen)
      if(nargin<2||isempty(name)||isempty(wid))
        error('must provide name and width'); 
      end
      if(nargin<3||isempty(rlen)), rlen = 1; end
      if(numel(wid)==1), wid = [1;2]*wid; end
      r = wid(2)/wid(1);
      subinfo.device = ['inv_',num2str(r)];
      subinfo.I_factor = [0;1];
      subinfo.C_factor = [1;1]*(1+r);
      % mode accept only one width for nmos
      this = this@coho_leaf(subinfo,name,wid(1),rlen);  
      this.i = this.add_port(node('i')); 
      this.o = this.add_port(node('o'));
      this.finalize;
    end
  end
end
