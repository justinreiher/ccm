% This circuit defines a Inverter-Ring-Oscillator circuit.  
%
% A 'IRO' has nodes 'xs', which is a N*1 cell. Each element is the output of the i-th stage 
%
% To create a circuit, requires
%   1. name: circuit name
%   2. N:    number of stage, must be odd number
%   3. wid:  circuit width
%            wid(1) is the inverter width
%   4. rlen: relative circuit length, use 1 by deault.
% 
% e.g. c = IRO('iro_5',5,1e-5);
%
classdef IRO < circuit
  properties (GetAccess='public', SetAccess='private');
    xs = {};
  end 
  methods
    function this = IRO(name,N,wid,rlen) 
      if(nargin<3||isempty(name)||isempty(N)||isempty(wid))
        error('not enough parameters'); 
      end
      if(nargin<4||isempty(rlen)), rlen = 1; end

      % check parameters
      if(mod(N,2)~=1)
        error('N must be an odd number');
      end

      this = this@circuit(name);

      % create nodes
      for i=1:N
        this.xs{i} = this.add_port(node(['x_',num2str(i)]));
      end

      % create basic elements
      invs = cell(N,1);
      for i=1:N
        invs{i} = inverter(['inv_',num2str(i)],wid,rlen);
        this.add_element(invs{i}); 
      end

      % connect elements
      for i=1:N-1
        this.connect(this.xs{i},invs{i}.o,invs{i+1}.i); 
      end
      this.connect(this.xs{N},invs{N}.o,invs{1}.i); 

      this.finalize;
    end

  end
end

