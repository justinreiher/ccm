% This class defines a circuit node
classdef node < matlab.mixin.Copyable %handle
  properties (GetAccess='public', SetAccess='private');
    name=[];
  end
  methods
    function this = node(name)
      this.name = name;
    end
  end
end
