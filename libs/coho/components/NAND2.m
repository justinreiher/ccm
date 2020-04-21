% This circuit defines a NAND circuit. 
% A 'NAND' has nodes: 
%   1. i1/i2: first and second inputs
%   2. o:     output
% To create a circuit, requires
%   1. name: circuit name
%   2. wid:  circuit width
%            wid(1) is the nmos width
%            wid(2) is the pmos width, use wid(1) by default
%   3. rlen: relative circuit length, use 1 by deault.
% E.g. n = NAND('nand',[1e-5;1e-5]);
classdef NAND2 < circuit
  properties (GetAccess='public', SetAccess='private');
    i1,i2; vdd; gnd; % input
    o; %output
  end 
  methods
    function this = NAND2 (name,wid,rlen) 
      if(nargin<2||isempty(name)||isempty(wid))
        error('must provide name and width'); 
      end
      if(nargin<3||isempty(rlen)), rlen = 1; end
      if(numel(wid)==1), wid=[1;1;1;1]*wid; end
      this = this@circuit(name); 
      
      assert(length(wid) == 4);
      
      widN1 = wid(1);
      widN2 = wid(2);
      widP1 = wid(3);
      widP2 = wid(4);
      
      this.vdd = this.add_port(node('vdd'));
      this.gnd = this.add_port(node('gnd'));
      
      this.i1  = this.add_port(node('i1'));
      this.i2  = this.add_port(node('i2'));
      this.o   = this.add_port(node('o'));
      x0        = this.add_port(node('x'));
      
      n1 = nmos('txN_1','wid',widN1,'rlen',rlen); this.add_element(n1);
      n2 = nmos('txN_2','wid',widN2,'rlen',rlen); this.add_element(n2);
      p1 = pmos('txP_1','wid',widP1,'rlen',rlen); this.add_element(p1);
      p2 = pmos('txP_2','wid',widP2,'rlen',rlen); this.add_element(p2);

      this.connect(this.vdd,p1.s,p2.s,p1.b,p2.b);
      this.connect(this.gnd,n2.s,n1.b,n2.b);
      this.connect(this.i1,n1.g,p1.g);
      this.connect(this.i2,n2.g,p2.g);
      this.connect(this.o,n1.d,p1.d,p2.d);
      this.connect(x0,n1.s,n2.d);
      this.finalize;
    end
  end
end
