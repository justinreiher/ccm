classdef rro_stage < circuit
  properties (GetAccess='public', SetAccess='private');
    i1,i2,o1,o2;
  end 
  methods
    % options
    %   wid: wid(1) is the size of forward inveter, wid(2) is the size of cross-couple INV 
    function this = rro_stage(name,wid,rlen) 
      if(nargin<2||isempty(name)||isempty(wid))
        error('not enough parameters'); 
      end
      if(nargin<3||isempty(rlen)), rlen = 1; end
      this = this@circuit(name);

      % create nodes
      this.i1 = this.add_port(node('i1')); this.i2 = this.add_port(node('i2'));
      this.o1 = this.add_port(node('o1')); this.o2 = this.add_port(node('o2'));

      % create elements
      fwdInv1 = INV('fwdInv1',wid(1),rlen); this.add_element(fwdInv1);
      fwdInv2 = INV('fwdInv2',wid(1),rlen); this.add_element(fwdInv2);
      ccInv   = rro_ccInv('ccInv',wid(2),rlen);  this.add_element(ccInv);

      % connect elements
      this.connect(this.i1,fwdInv1.i); this.connect(this.i2,fwdInv2.i);
      this.connect(this.o1,fwdInv1.o,ccInv.x1); this.connect(this.o2,fwdInv2.o,ccInv.x2); 
      this.finalize;
    end
  end
end

