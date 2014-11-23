% This class defines n-stage rambus ring oscillator
% use to test flatten.
classdef rambusOsc < circuit
  properties (GetAccess='public', SetAccess='private');
    xs = {};
  end 
  methods
    % options
    %   'N': number of stages
    %   wid: wid(1) is the size of forward inveter, wid(2) is the size of cross-couple inverter 
    function this = rambusOsc(name,N,wid,rlen) 
      if(nargin<3||isempty(name)||isempty(N)||isempty(wid))
        error('not enough parameters'); 
      end
      if(nargin<4||isempty(rlen)), rlen = 1; end

      % check parameters
      if(mod(N,2)~=0)
        error('N must be an odd number');
      end

      this = this@circuit(name);

      % create nodes
      for i=1:N
        this.xs{i,1} = this.add_port(node(['x_u',num2str(i)]));
        this.xs{i,2} = this.add_port(node(['x_b',num2str(i)]));
      end
      % create basic elements
      stages = cell(N); 
      for i=1:N
        stages{i} = rro_stage(['s_',num2str(i)],wid,rlen); 
        this.add_element(stages{i});
      end
      % connect elements
      for i=1:N-1
        this.connect(this.xs{i,1},stages{i}.o1,stages{i+1}.i1);
        this.connect(this.xs{i,2},stages{i}.o2,stages{i+1}.i2);
      end
      % swap at the end
      this.connect(this.xs{N,1},stages{N}.o1,stages{1}.i2);
      this.connect(this.xs{N,2},stages{N}.o2,stages{1}.i1);
      this.finalize;
    end
  end
end

