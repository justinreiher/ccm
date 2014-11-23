classdef rro_ccInv < circuit
  properties (GetAccess='public', SetAccess='private');
    x1,x2
  end 
  methods
    % options
    %   wid: wid(1) is the size of forward inveter, wid(2) is the size of cross-couple INV 
    function this = rro_ccInv(name,wid,rlen) 
      if(nargin<2||isempty(name)||isempty(wid))
        error('not enough parameters'); 
      end
      if(nargin<3||isempty(rlen)), rlen = 1; end
      this = this@circuit(name);

      % create nodes
      this.x1 = this.add_port(node('x1')); this.x2 = this.add_port(node('x2')); 

      % create elements
      inv1 = INV('inv1',wid,rlen); this.add_element(inv1);
      inv2 = INV('inv2',wid,rlen); this.add_element(inv2);

      % connect elements
      this.connect(this.x1,inv1.i,inv2.o); 
      this.connect(this.x2,inv2.i,inv1.o); 
      this.finalize;
    end
  end
end
