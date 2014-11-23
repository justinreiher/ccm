classdef lib_coho_test_circuit_isrc < circuit
  properties (GetAccess='public', SetAccess='private');
    x;
  end 
  methods
    function this = lib_coho_test_circuit_isrc (name,i,c) 
      if(nargin<3), error('not enough parameter'); end
      this = this@circuit(name); 
      this.x = this.add_port(node('x'));
       
      src = isrc('src',i,true); this.add_element(src); 
      cap = capacitor('cap',c,true); this.add_element(cap); 
      
      this.connect(this.x,src.p,cap.x); 
      this.finalize;
    end
  end
end
