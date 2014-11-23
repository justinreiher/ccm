% This circuit defines a weak feedback celement circuit. 
% Detail information available on the paper 'Formal Verification of C-element circuits' 
% A 'celement' has three nodes: 
%   1. i1/i2: input
%   2. o: output 
% To create a circuit, requires
%   1. name: circuit name
%   2. wid:  circuit width 
%            wid(1) is the width for nmos 
%            wid(2) is the width for pmos, use 2*wid(1) by default 
%            wid(3) is the width for output inverter, use 0.5*wid(1) by default
%            wid(4) is the width for keeper inverter, use 0.5*wid(3) by default
%   3. rlen: relative circuit length
% E.g. c = CELEMENT('celement',[1;2;0.5;0.25]*1e-5,1)
classdef CELEMENT < circuit
  properties (GetAccess='public', SetAccess='private');
    i1; i2; % input
    o; % output
  end 
  methods
    function this = CELEMENT(name,wid,rlen) 
      if(nargin<2||isempty(name)||isempty(wid))
        error('must provide name and width'); 
      end
      if(nargin<3||isempty(rlen)), rlen = 1; end
      if(numel(wid)<2), wid = [1;2]*wid; end
      if(numel(wid)<3), wid(3) = 0.5*wid(1); end
      if(numel(wid)<4), wid(4) = 0.5*wid(3); end

      this = this@circuit(name); 

      % nodes
      this.i1 = this.add_port(node('i1'));
      this.i2 = this.add_port(node('i2'));
      this.o  = this.add_port(node('o'));
      xn      = this.add_port(node('xn'));
      xp      = this.add_port(node('xp'));
      xy      = this.add_port(node('xy'));

      % elements
      n1 = nmos('n1',wid(1),1,rlen); this.add_element(n1);
      p1 = pmos('p1',wid(2),1,rlen); this.add_element(p1);
      n2 = nmos('n2',wid(1),0,rlen); this.add_element(n2);
      p2 = pmos('p2',wid(2),0,rlen); this.add_element(p2);
      oinv = inverter('oinv',wid(3),rlen); this.add_element(oinv);
      kinv = inverter('kinv',wid(4),rlen); this.add_element(kinv);

      % connection
      this.connect(this.i1, n1.g, p1.g); 
      this.connect(this.i2, n2.g, p2.g); 
      this.connect(this.o, oinv.o, kinv.i); 
      this.connect(xn, n1.d, n2.s); 
      this.connect(xp, p1.d, p2.s); 
      this.connect(xy, n2.d, p2.d, oinv.i, kinv.o); 
     
      this.finalize;
    end
  end
end
