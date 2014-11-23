% This circuit defines a D_LATCH circuit using NANDs. 
% A 'latch' has three nodes: 
%   1. d: input
%   2. q: output
%   3. en: enable
% To create a circuit, requires
%   1. name: circuit name
%   2. wid:  circuit width 
%            wid(1) is the width for NAND
%            wid(2) is the width for the inverter, use wid(1)/2 by default 
%   3. rlen: relative circuit length, use 1 by deault.
% E.g. c = LATCH('latch',[1e-5;0.5e-5]);
classdef LATCH < circuit
  properties (GetAccess='public', SetAccess='private');
    d,q,en;
  end 
  methods
    function this = LATCH (name,wid,rlen) 
      if(nargin<2||isempty(name)||isempty(wid))
        error('must provide name and width'); 
      end
      if(nargin<3||isempty(rlen)), rlen = 1; end
      if(numel(wid)==1), wid=[1;0.5]*wid; end
      this = this@circuit(name); 
      this.d = this.add_port(node('d')); 
      this.q = this.add_port(node('q')); 
      this.en = this.add_port(node('en'));
      dbar   = this.add_port(node('dbar'));
      qbar   = this.add_port(node('qbar'));
      x1     = this.add_port(node('x1'));
      x2     = this.add_port(node('x2'));
      
      nand_s0 = NAND('NAND_S0',wid(1),rlen); this.add_element(nand_s0); 
      nand_r0 = NAND('NAND_R0',wid(1),rlen); this.add_element(nand_r0); 
      nand_s = NAND('NAND_S',wid(1),rlen); this.add_element(nand_s); 
      nand_r = NAND('NAND_R',wid(1),rlen); this.add_element(nand_r); 
      inv    = inverter('inv',wid(2),rlen); this.add_element(inv);

      this.connect(this.d, inv.i, nand_s0.i1); 
      this.connect(dbar, inv.o, nand_r0.i1);
      this.connect(this.en, nand_s0.i2, nand_r0.i2);
      this.connect(x1, nand_s0.o, nand_s.i1);
      this.connect(x2, nand_r0.o, nand_r.i1);
      this.connect(this.q, nand_s.o, nand_r.i2);
      this.connect(qbar, nand_r.o, nand_s.i2);
      this.finalize;
    end
  end
end
