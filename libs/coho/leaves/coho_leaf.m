% This class defines an interface for basic elements used in COHO
% For a COHO leaf circuit, the constructor usually requires
%   1. name: circuit name, 
%   2. wid:  circuit width
%   3. rlen: relative circuit length (compared to the minimum length)
% This class provides all functions to support all three interfaces of 
% the 'circuit' abstract class. The leaf circuit class only need to 
% implement the constructor, which needs to pass the 'subinfo' with
%   1. 'device': device name corresponds to the 'mat' file
%   2. I_factor: how I is calculated for all nodes
%   2. C_factor: how C is calcualted for all nodes

% TODO. We do not support multiple models now. 
% It is easy to implement. However, the interface of the function 
% must be changed (array -> cell). I hate this. Update it later.
% TODO: remove vars from I_ldi or optimize for it?
classdef coho_leaf< circuit
  properties (Constant)
    capPerUnit = 2e-9; 
  end
  properties (GetAccess='public', SetAccess='private');
    subinfo; 
    params; 
  end 

  methods
    % device:  corresponding mat file name 
    % name:    name of the circuit
    % wid:     wid of the device
    % rlen:    relative length compare with the mimimum length 
    function this = coho_leaf(subinfo,name,params)
      if(nargin<3), error('must provide subinfo, name and parameters'); end
      this = this@circuit(name);
      this.subinfo.device= subinfo.device; 
      this.subinfo.I_factor = subinfo.I_factor(:); % make it col vector
      this.subinfo.C_factor = subinfo.C_factor(:);
      this.params = params; 
    end

    % support simulation
    function is = ifc_simu(this)
      is = true;
    end

    % v is NxM, i is NxM
    function i = I(this,v,varargin)
      i = interp_device(this.subinfo.device,v,this.wid,this.rlen);  % col vector
      i = repmat(i',this.nodeNum,1); % each col for all nodes
      factor = repmat(this.subinfo.I_factor,1,size(i,2));
      i = i.*factor;  
    end

    % NOTE: assume cap does not depends on voltage
    % NOTE: assume cap is linear with wid, independent of length? 
    function c = C(this,v,varargin) 
      if(nargin<2||isempty(v)), v = zeros(this.nodeNum,1); end
      c = this.capPerUnit*this.params.wid*this.params.rlen; %% TODO: check if it's correct? 
      c = repmat(c,this.nodeNum,size(v,2)); % each col for all nodes
      factor = repmat(this.subinfo.C_factor,1,size(c,2));
      c = c.*factor; 
    end

    % support small signal analysis 
    function is = ifc_ssan(this)
      is = true;
    end

    % v is Vx1; didv is NxN 
    function didv = dIdV(this,v,varargin) 
      assert(size(v,2)==1)
      didv = interp_device_jac(this.subinfo.device,v,this.param.wid,this.param.rlen); 
      didv = repmat(didv',this.nodeNum,1); % didv(i,j) = dI_i/dV_j 
      factor = repmat(this.subinfo.I_factor,1,size(didv,2));
      didv = didv.*factor;
    end

    % support verification
    function is = ifc_verify(this)
      is = true;
    end

    % Region: an object of coho_shape  class
    function [c,err] = I_ldi(this,region,vars,varargin)
      if(nargin<3||isempty(vars)), vars = (1:this.nodeNum); end
      [c,err] = lft_device(this.subinfo.device,region,this.wid,this.rlen); 
      c = repmat(c',this.nodeNum,1); err = repmat(err,this.nodeNum,1); % repeat for all nodes
      % c is Nx(N+1), err is Nx1
      factor = this.subinfo.I_factor; % Nx1
      c = c.*repmat(factor,1,size(c,2)); err = err.*abs(factor);
      % select variables
      c = c(vars,:); err = err(vars);
    end
  end % method

  methods(Access=protected)
    function id = vectorize_subtype_id(this)
      id = this.subinfo.device;
    end
  end % methods

  methods(Static)
    function I = I_objs(Objs,V,varargin)
      % mrg: N = number of terminal per device;
      %      P = number of phase-space points at which to evaluate the current;
      %      K = number of devices with this device model
      [N,P,K] = size(V); assert(length(Objs)==K);
      device = Objs{1}.subinfo.device; I_factor = Objs{1}.subinfo.I_factor;
      % collect all width/rlength
      wid = zeros(K,1);  rlen = zeros(K,1); 
      for ind=1:K
        wid(ind) = Objs{ind}.wid; rlen(ind) = Objs{ind}.rlen;
      end
      % make V from N*P*K to N*(P*K)
      V = reshape(V,N,[]);
      % wid/rlen should b (P*K) x 1, same for all P
      wid = reshape(repmat(wid,1,P)',P*K,1); % [repmat(w1,P) ... repmat(wk,P)]
      rlen = reshape(repmat(rlen,1,P)',P*K,1); 
      % call vectorized I function 
      i = interp_device(device,V,wid,rlen);   % calculate i 1x(P*K)
      i = repmat(i',N,1);                  % repeat for all nodes: N*(P*K)
      factor = repmat(I_factor,1,P*K);     % select
      i = i.*factor;                       % N*(P*K)
      % change size back to N*P*K
      I = reshape(i,[N,P,K]);
    end

    function C = C_objs(Objs,V,varargin)
      [N,P,K] = size(V); assert(length(Objs)==K); 
      device = Objs{1}.subinfo.device; C_factor = Objs{1}.subinfo.C_factor;
      % collect all width/rlength
      wid = zeros(K,1);  rlen = zeros(K,1); 
      for i=1:K
        wid(i) = Objs{i}.wid; rlen(i) = Objs{i}.rlen;
      end
      % make V from N*P*K to N*(P*K)
      V = reshape(V,N,[]);
      % wid/rlen should b (P*K) x 1
      wid = reshape(repmat(wid,1,P)',P*K,1); % [repmat(w1,P) ... repmat(wk,P)]
      rlen = reshape(repmat(rlen,1,P)',P*K,1); 
      % call vectorized cap function 
      cap = Objs{1}.capPerUnit.*wid.*rlen;            % calculate cap 
      cap = repmat(cap',N,1);                % repeat for all nodes 
      factor = repmat(C_factor,1,P*K);     % select
      cap = cap.*factor;                     % N*(P*K)
      % change back to N*P*K
      C = reshape(cap,[N,P,K]);
    end
  end % method
end
