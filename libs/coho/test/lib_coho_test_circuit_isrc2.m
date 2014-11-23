classdef lib_coho_test_circuit_isrc2 < circuit
  properties (GetAccess='public', SetAccess='private');
    x,y;
  end 
  methods
    function this = lib_coho_test_circuit_isrc2 (name,i,c) 
      if(nargin<3), error('not enough parameter'); end
      this = this@circuit(name); 
      this.x = this.add_port(node('x'));
      this.y = this.add_port(node('y'));
       
      src = isrc('src',i,false); this.add_element(src); 
      cap = capacitor('cap',c,false); this.add_element(cap); 
      
      this.connect(this.x,src.p,cap.x); 
      this.connect(this.y,src.m,cap.y); 
      this.finalize;
    end
  end
end
