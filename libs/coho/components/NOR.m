% This circuit defines a NOR circuit. 
% A 'NOR' has nodes: 
%   1. i1/i2: first and second inputs
%   2. o:     output
% To create a circuit, requires
%   1. name: circuit name
%   2. wid:  circuit width
%            wid(1) is the nmos width
%            wid(2) is the pmos width, use 4*wid(1) by default
%   3. rlen: relative circuit length, use 1 by deault.
% E.g. n = NOR('nor',[1e-5;4e-5]);
classdef NOR < circuit
  properties (GetAccess='public', SetAccess='private');
    i1,i2; o; % input/output
  end 
  properties (Constant)
    cap_factor = -0.75;
  end
  methods
    function this = NOR (name,wid,rlen) 
      if(nargin<2||isempty(name)||isempty(wid))
        error('must provide name and width'); 
      end
      if(nargin<3||isempty(rlen)), rlen = 1; end
      if(numel(wid)==1), wid=[1;4]*wid; end
      this = this@circuit(name);
      this.i1  = this.add_port(node('i1'));
      this.i2  = this.add_port(node('i2'));
      this.o   = this.add_port(node('o'));
      x        = this.add_port(node('x'));
      
      n1 = nmos('n1',wid(1),1,rlen); this.add_element(n1);
      n2 = nmos('n2',wid(1),1,rlen); this.add_element(n2);
      p1 = pmos('p1',wid(2),1,rlen); this.add_element(p1);
      p2 = pmos('p2',wid(2),0,rlen); this.add_element(p2);
      % TODO modify cap function?
      c1 = p1.C; c2 = p2.C;
      cc = (c1(p1.find_port_index(p1.g))+c2(p2.find_port_index(p2.d))); % source and gain
      c = capacitor('c',this.cap_factor*cc,1); this.add_element(c);

      this.connect(this.i1,n1.g,p1.g);
      this.connect(this.i2,n2.g,p2.g);
      this.connect(this.o,n1.d,n2.d,p2.d);
      this.connect(x,p1.d,p2.s,c.x); 
      this.finalize;
    end
  end
end
