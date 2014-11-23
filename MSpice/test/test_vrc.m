classdef test_vrc < circuit
  properties (GetAccess='public', SetAccess='private');
    x,y,z; 
  end 
  methods
    function this = test_vrc (name,v,r,c) 
      if(nargin<4), error('not enough parameter'); end
      this = this@circuit(name); 
      this.x = this.add_port(node('x'));
      this.y = this.add_port(node('y'));
      this.z = this.add_port(node('z'));
       
      src = vsrc('vsrc',v,false); this.add_element(src); 
      res = resistor('resistor',r,false); this.add_element(res); 
      cap = capacitor('cap',c,false); this.add_element(cap); 
      
      this.connect(this.x,src.p,res.x); 
      this.connect(this.y,res.y,cap.x); 
      this.connect(this.z,src.m,cap.y);

      this.finalize;
    end
  end
end
