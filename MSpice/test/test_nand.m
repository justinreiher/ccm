classdef test_nand < circuit
  properties (GetAccess='public', SetAccess='private');
    i1,i2; o; % input/output
    vdd,gnd;
  end 
  properties (Constant)
    cap_factor = -0.75;
  end
  methods
    function this = test_nand (name,wid,rlen) 
      if(nargin<2||isempty(name)||isempty(wid))
        error('must provide name and width'); 
      end
      if(nargin<3||isempty(rlen)), rlen = 1; end
      if(numel(wid)==1), wid=[1;1]*wid; end
      this = this@circuit(name); 
      this.vdd = this.add_port(node('vdd'));
      this.gnd = this.add_port(node('gnd'));
      this.i1  = this.add_port(node('i1'));
      this.i2  = this.add_port(node('i2'));
      this.o   = this.add_port(node('o'));
      x        = this.add_port(node('x'));
      
      n1 =  nmos('n1',wid(1),false,rlen); this.add_element(n1);
      n2 =  nmos('n2',wid(1),false,rlen); this.add_element(n2);
      p1 =  pmos('p1',wid(2),false,rlen); this.add_element(p1);
      p2 =  pmos('p2',wid(2),false,rlen); this.add_element(p2);

      %  reduce capacitance of internal nodes
      c1 = n1.C; c2 = n2.C;
      cc = (c1(n1.find_port_index(n1.g))+c2(n2.find_port_index(n2.d))); % source and gain
      c = capacitor('c',this.cap_factor*cc,1); this.add_element(c);
      
      this.connect(this.vdd,p1.s,p2.s);
      this.connect(this.gnd,n2.s); 
      this.connect(this.i1,n1.g,p1.g);
      this.connect(this.i2,n2.g,p2.g);
      this.connect(this.o,n1.d,p1.d,p2.d);
      this.connect(x,n1.s,n2.d,c.x); 
      this.finalize;
    end
  end
end
