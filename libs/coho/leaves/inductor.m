% Inductor will make the ODE second-order, so do not support now
classdef inductor < circuit 
  methods
    function this = inductor(name,l)
      error('we do not support inductor now');
    end
  end
end
