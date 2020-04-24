% Definition of a THREE_RING_OSC a three stage ring oscillator with
% INV4 elements
%
% The three ring oscillator has 2 input nodes:
% 1. vdd:    power supply source
% 2. gnd:    circuit ground
% The three ring oscillator has 3 output nodes:
% 1. v1:     the first output of the oscillator
% 2. v2:     the second output of the oscillator
% 3. v3:     the third output of the oscillator
%
% To create a THREE_RING_OSC, requires
% 1. name: oscillator name
% 2. wid: circuit width
%       - if one wid is given that width is applied to all devices
%       - otherwise width needs to be of length 6:
%           wid(1:2) = nmos inv 1, pmos inv 1
%           wid(3:4) = nmos inv 2, pmos inv 2
%           wid(5:6) = nmos inv 3, pmos inv 3
% 3. rlen: relative device length, use 1 by default.
% E.g. myOsc = THREE_RING_OSC('osc0',450e-7,1)
classdef THREE_RING_OSC < circuit
    
    properties (GetAccess = 'public', SetAccess = 'private')
        vdd,gnd; %Input
    end
    
    methods
        function this = THREE_RING_OSC(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1)
                widInv0=[1;1]*wid;
                widInv1=[1;1]*wid;
                widInv2=[1;1]*wid;
            else
                assert(length(wid) == 6)
                widInv0 = wid(1:2);
                widInv1 = wid(3:4);
                widInv2  = wid(5:end);
            end
            this = this@circuit(name);
            %set the relative lengths by default to be be the same
            rlen0 = rlen;
            rlen1 = rlen;
            rlen2 = rlen;
            if(iscell(rlen))
                rlen0 = rlen{1};
                rlen1 = rlen{2};
                rlen2 = rlen{3};
            end
            %define circuit nodes
            %inputs
            this.vdd = this.add_port(node('vdd'));
            this.gnd = this.add_port(node('gnd'));
            
            %internal node/outputs
            v0 = this.add_port(node('v0'));
            v1 = this.add_port(node('v1'));
            v2 = this.add_port(node('v2'));
            
            inv0 = INV4('inv0',widInv0,rlen0); this.add_element(inv0);
            inv1 = INV4('inv1',widInv1,rlen1); this.add_element(inv1);
            inv2 = INV4('inv2',widInv2,rlen2); this.add_element(inv2);
            
            this.connect(v2,inv2.o,inv0.i);
            this.connect(v0,inv0.o,inv1.i);
            this.connect(v1,inv1.o,inv2.i);
            this.connect(this.vdd,inv0.vdd,inv1.vdd,inv2.vdd);
            this.connect(this.gnd,inv0.gnd,inv1.gnd,inv2.gnd);
            
            
            this.finalize;
            
        end
        
    end
end

