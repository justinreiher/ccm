% This circuit defines a Rambus-Ring-Oscillator circuit. 
%
% A 'RRO' has nodes 'xs', which is a N*2 cell 
%   1. xs{i,1}: the up nodes for the i-th stage. 
%   2. xs{i,2}: the bottom nodes for the i-th stage. 
%
% To create a circuit, requires
%   1. name: circuit name
%   2. N:    number of stage, must be even number
%   3. wid:  circuit width
%            wid(1) is the forward inverter width
%            wid(2) is the cross-coupled inverter width, use wid(1) by default
%   4. rlen: relative circuit length, use 1 by deault.
% 
% e.g. c = RRO('rro_4',4,[1;1]*1e-5);
%
classdef RRO < circuit
  properties (GetAccess='public', SetAccess='private');
    xs = {};
  end 
  methods
    function this = RRO(name,N,wid,rlen) 
      if(nargin<3||isempty(name)||isempty(N)||isempty(wid))
        error('not enough parameters'); 
      end
      if(nargin<4||isempty(rlen)), rlen = 1; end
      if(numel(wid)==1), wid=[1;1]*wid; end

      % check parameters
      if(mod(N,2)~=0)
        error('N must be an even number');
      end

      this = this@circuit(name);

      % create nodes
      for i=1:N
        this.xs{i,1} = this.add_port(node(['x_u_',num2str(i)]));
        this.xs{i,2} = this.add_port(node(['x_b_',num2str(i)]));
      end

      % create basic elements
      fwd_inv = cell(N,2); cc_inv = cell(N,2);
      for i=1:N
        fwd_inv{i,1} = inverter(['fwd_u_inv_',num2str(i)],wid(1),rlen); 
        fwd_inv{i,2} = inverter(['fwd_b_inv_',num2str(i)],wid(1),rlen); 
        cc_inv{i,1} = inverter(['cc_u_inv_',num2str(i)],wid(2),rlen); 
        cc_inv{i,2} = inverter(['cc_b_inv_',num2str(i)],wid(2),rlen); 
        this.add_element(fwd_inv{i,1}); this.add_element(fwd_inv{i,2});
        this.add_element(cc_inv{i,1}); this.add_element(cc_inv{i,2});
      end

      % connect elements
      for i=1:N-1
        this.connect(this.xs{i,1},fwd_inv{i,1}.o,fwd_inv{i+1,1}.i, cc_inv{i,1}.o, cc_inv{i,2}.i); 
        this.connect(this.xs{i,2},fwd_inv{i,2}.o,fwd_inv{i+1,2}.i, cc_inv{i,2}.o, cc_inv{i,1}.i); 
      end
      % swap at the end
      this.connect(this.xs{N,1},fwd_inv{N,1}.o,fwd_inv{1,2}.i, cc_inv{N,1}.o, cc_inv{N,2}.i); 
      this.connect(this.xs{N,2},fwd_inv{N,2}.o,fwd_inv{1,1}.i, cc_inv{N,2}.o, cc_inv{N,1}.i); 

      this.finalize;
    end

  end
end

