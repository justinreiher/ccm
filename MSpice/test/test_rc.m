classdef test_rc < circuit
  properties (GetAccess='public', SetAccess='private');
    x;
  end 
  methods
    function this = test_rc (name,r,c) 
      if(nargin<3), error('not enough parameter'); end
      this = this@circuit(name); 
      this.x = this.add_port(node('x'));
       
      res = resistor('resistor',r,true); this.add_element(res); 
      cap = capacitor('cap',c,true); this.add_element(cap); 
      
      this.connect(this.x,res.x,cap.x); 

      this.finalize;
    end
  end
end
