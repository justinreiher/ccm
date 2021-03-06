% Definition of a TWO_FLOP_STRONG_ARM_SYNC with STRONG_ARM_FF and INV4 elements
% defines a two flip-flop strong arm synchronizer
% The strong arm synchronizer has 4 input nodes:
% 1. vdd:    power supply source
% 2. gnd:    circuit ground
% 3. d:      data input
% 4. clk:    clock input
%
% The Strong Arm synchronizer has 2 output nodes:
% 1. q:      output
% 2. qbar:   output bar (i.e. the opposite polarity to q)
%
% To create a TWO_FLOP_STRONG_ARM_SYNC, requires
% 1. name: strong arm synchronizer name
% 2. wid: circuit width, wid(1:2)   is the input buffer for the synchronizer
%                                   (INV4)
%                        wid(3:4)   is the logical inverse of the buffered
%                                   input (INV4)
%                        wid(5:22)  is for the Master stage of the 
%                                   synchronizer (STRONG_ARM_FF)
%                        wid(23:40) is for the Slave stage of the 
%                                   synchronizer (STRONG_ARM_FF)
%
% 3. rlen: relative circuit length, use 1 by default.
% E.g. sync = TWO_FLOP_STRONG_ARM_SYNC('saSync',450e-7,1)

classdef TWO_FLOP_STRONG_ARM_SYNC < circuit
    
    properties (GetAccess = 'public', SetAccess = 'private')
        vdd,gnd,d,clk; q,qbar; %Input ; Output
    end
    
    methods
        function this = TWO_FLOP_STRONG_ARM_SYNC(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1)
                widInv0=[1;1]*wid;
                widInv1 = [1;1]*wid;
                widMaster = wid;
                widSlave = wid;
            else
                assert(length(wid) == 40)
                widInv0 = wid(1:2);
                widInv1 = wid(3:4);
                widMaster = wid(5:22);
                widSlave  = wid(23:end);

            end
            this = this@circuit(name);
            %set the relative lengths by default to be be the same
            rlenInv0 = rlen;
            rlenInv1 = rlen;
            rlenMaster = rlen;
            rlenSlave = rlen;
            
            if(iscell(rlen))
                rlenInv0 = rlen{1};
                rlenInv1 = rlen{2};
                rlenMaster = rlen{3};
                rlenSlave = rlen{4};
            end
            %define circuit nodes
            %inputs
            this.d = this.add_port(node('d'));
            this.clk = this.add_port(node('clk'));
            this.vdd = this.add_port(node('vdd'));
            this.gnd = this.add_port(node('gnd'));
            
            %internal node
            dinBar = this.add_port(node('dinBar'));
            din    = this.add_port(node('din'));
            %outpus
            this.q = this.add_port(node('q'));
            this.qbar = this.add_port(node('qbar'));
            
            inv0     = INV4('bufferDin',widInv0,rlenInv0); this.add_element(inv0);
            inv1     = INV4('bufferDinBar',widInv1,rlenInv1); this.add_element(inv1);
            sARM_ffMaster = STRONG_ARM_FF('saFF_0',widMaster,rlenMaster); this.add_element(sARM_ffMaster);
            sARM_ffSlave = STRONG_ARM_FF('saFF_1',widSlave,rlenSlave); this.add_element(sARM_ffSlave);

            
            this.connect(this.d,inv0.i);
            this.connect(dinBar,inv0.o,inv1.i,sARM_ffMaster.dbar);
            this.connect(din,inv1.o,sARM_ffMaster.d);
            this.connect(sARM_ffMaster.q,sARM_ffSlave.d);
            this.connect(sARM_ffMaster.qbar,sARM_ffSlave.dbar);
            this.connect(this.clk,sARM_ffMaster.clk,sARM_ffSlave.clk);
            this.connect(this.vdd,sARM_ffMaster.vdd,sARM_ffSlave.vdd,inv0.vdd,inv1.vdd);
            this.connect(this.gnd,sARM_ffMaster.gnd,sARM_ffSlave.gnd,inv0.gnd,inv1.gnd);
            
            this.connect(this.q,sARM_ffSlave.q);
            this.connect(this.qbar,sARM_ffSlave.qbar);
            this.finalize;
            
        end
        
    end
end

