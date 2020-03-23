% This class defines the 'inverter' circuit for COHO library 
% An 'inverter' has two nodes: i(input) and o(output)
% To create an inveter, need to provide
%   1. name: circuit name
%   2. wid:  circuit width 
%            wid(1) is the nmos width
%            wid(2) is the pmos width, use 2*wid(1) if not provided.
%   3. rlen: (circuit length)/(minimum length), use 1 by default
% E.g. inv = INV('inv',[1e-5,2e-5],1)
% NOTE: It's different from 'inverter' that 'INV' consists of 'nmos' and 'pmos'. 
%  i.e. a function call of 'INV' requires two function calls from 'nmos','pmos'. 
%  while a function call of 'inverter' issues one call from the mat file directly. 
classdef INV4 < circuit
  properties (GetAccess='public', SetAccess='private');
    vdd,gnd,i; o;
  end 
  methods
    % wid: wid(1) is the width of nmos, and wid(2) is the width of pmos
    %      if wid(2) is not provided, 2*wid(1) is used by default
    % rlen: relative transitor length, must be scaler 
    function this = INV4 (name,wid,rlen)
      if(nargin<2||isempty(name)||isempty(wid))
        error('must provide name and width'); 
      end
      if(nargin<3||isempty(rlen)), rlen = 1; end
      if(numel(wid)==1), wid = [1;1]*wid; end
      this    = this@circuit(name);
      
      this.vdd = this.add_port(node('vdd'));
      this.gnd = this.add_port(node('gnd'));
      this.i  = this.add_port(node('i'));
      this.o  = this.add_port(node('o'));
      
      %create devices wid is of the form wid = [widN,widP] if wid is
      %called as a single number then both the P and N device will
      %have the same width.
      
      widN = wid(1);
      widP = wid(2);
      
      n = nmos('n','wid',widN,'rlen',rlen); this.add_element(n);
      p = pmos('p','wid',widP,'rlen',rlen); this.add_element(p); 

      this.connect(this.i,n.g,p.g);
      this.connect(this.o,n.d,p.d);
      this.connect(this.vdd,p.s,p.b);
      this.connect(this.gnd,n.s,n.b);
      this.finalize;
    end
  end
end
