% This class defines the 'PMOS' for COHO library
% The circuit has four nodes: s,g,d,b
% To create a PMOS, need to provide
%   1. name: circuit name
%   2. wid:  circuit width
%   2. ctp:  source connected to vdd, false by default
%   3. rlen: (circuit length)/(minimum length), use 1 by default
% E.g. p = pmos('pmos',1e-5,0,true)
classdef pmos < coho_leaf
  properties (GetAccess='public', SetAccess='private'); 
    s; d; g; b; model;
  end
  methods
    function this = pmos(name,varargin)
      
      params = [];
      
      %figure out if the input is all numeric
      try
          testNumericInput = isnumeric([varargin{:}]);
      catch
          testNumericInput = false;
      end
       
       
      %get the number of arguments passed into the pmos function 
      [~,numArg] = size(varargin);
      if(numArg == 0)
          %if the number of arguments in varargin is 0 then insufficient
          %data was provided and an error is thrown
          error('pmos called with insufficient information')
      elseif(numArg == 1 && isstruct(varargin{1}))
          %if the number of arguments in varargin is 1 and it's argument is
          %a struct then that is taken as the parameters for the device
          params = varargin{1};
      elseif(testNumericInput)
          %if all the values of varargin are numeric, then it is assumed
          %that it follows the old calling convention.
          warning('pmos created using deprecated function call in circuit definition')
          %old convention call has been made which follows:
          %pmos(name,wid,ctp,rlen)
          %This call is now packaged into params and calls the
          %transistorFactory to create the appropriate transistor
          params.wid       = varargin{1};
          try
              ctgFlag = varargin{2};
          catch
              ctgFlag = false;
          end
          params.gndFlag    = ctgFlag; %this is a legacy flag
          params.rlen       = varargin{3};
      elseif(mod(numArg,2)~=0)
          %if the number of arguments in varargin is not even then there is
          %a (key,value) pair imbalance.
          error('pmos called with imbalanced (key,value) pairs')
      else
          %otherwise it is assumed that the call is of the right form and
          %the parameter structure is constructed by walking through the
          %list of varargin, if varargin is even but does not have follow
          %the (string,value) pair convention then the default MATLAB error
          %is thrown. MATLAB complains about the dynamic nature of setfield
          %in this way, but this only happens at circuit construction.
          for i =1:(numArg/2)
              index = 2*i-1;
              key = varargin{index};
              value = varargin{index+1};
              params.(key) = value;
          end
      end
      
      %try to access a model from the collected parameters, if none exists
      %then assign the default model to be used for this device.
      try
          params.model;
      catch
          params.model = 'default';
      end
      
      subinfo.device='pmos';
      subinfo.I_factor=[1;0;-1;0];
      subinfo.C_factor=[1;1;1;0];    
      
      this = this@coho_leaf(subinfo,name,params);
      %transistor has 4 ports, drain, gate, source, body
      this.model = params.model;
      this.d = this.add_port(node('d'));
      this.g = this.add_port(node('g'));
      this.s = this.add_port(node('s'));
      this.b = this.add_port(node('b'));
      
      this.finalize;
    end
  end

  methods(Access=protected)
    function id = vectorize_subtype_id(this) 
        id = this.model;
    end
  end
end
