% This class define a shape class for ids_ldi function used by verification
classdef shape < matlab.mixin.Copyable %handle
  methods(Abstract)
    s = project(this,index)
    b = rec(this);
  end
end
