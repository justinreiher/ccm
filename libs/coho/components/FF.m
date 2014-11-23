% This circuit defines a D_FLIP_FLOP circuit using D_LATCH. 
% A 'flip-flop' is positive-edge-triggered, with three nodes: 
%   1. d: input
%   2. q: output
%   3. clk: clock 
% To create a circuit, requires
%   1. name: circuit name
%   2. wid:  circuit width 
%            wid(1) is the width for D_LATCH 
%            wid(2) is the width for the inverter, use wid(1)/2 by default 
%   3. rlen: relative circuit length, use 1 by deault.
% E.g. c = FF('flip-flop',[1e-5;0.5e-5]);

classdef FF < circuit
  properties (GetAccess='public', SetAccess='private');
    d,q,clk; 
  end 
  methods
    function this = FF(name,wid,rlen) 
      if(nargin<2||isempty(name)||isempty(wid))
        error('must provide name and width'); 
      end
      if(nargin<3||isempty(rlen)), rlen = 1; end
      if(numel(wid)==1), wid=[1;0.5]*wid; end

      this = this@circuit(name); 
      this.d = this.add_port(node('d')); 
      this.q = this.add_port(node('q')); 
      this.clk = this.add_port(node('clk'));
      clkbar = this.add_port(node('clkbar')); 
      x      = this.add_port(node('x'));
      
      latch_m = LATCH('LATCH_MASTER',wid(1),rlen); this.add_element(latch_m); 
      latch_s = LATCH('LATCH_SLAVE',wid(1),rlen); this.add_element(latch_s); 
      inv    = inverter('inv',wid(2),rlen); this.add_element(inv);

      this.connect(this.d, latch_m.d); 
      this.connect(this.clk, inv.i, latch_s.en);
      this.connect(clkbar, inv.o, latch_m.en);
      this.connect(x, latch_m.q, latch_s.d);
      this.connect(this.q, latch_s.q);

      this.finalize;
    end
  end
end
