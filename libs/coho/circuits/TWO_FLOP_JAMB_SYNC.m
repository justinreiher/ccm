% Definition of a TWO_FLOP_JAMB_SYNC with JAMB_FF and INV4 elements
%and defines a two flip-flop jamb latch synchronizers 
% The Jamb synchronizer has 6 input nodes:
% 1. vdd:    power supply source
% 2. gnd:    circuit ground
% 3. d:      data input
% 4. clk:    clock input
% 5. clkbar: clock bar (i.e. the opposite polarity to clock)
% 6. reset:  reset input to reset the jamb latches in the synchronizer
% The passgate latch has 2 output nodes:
% 1. q:      output
% 2. qbar:   output bar (i.e. the opposite polarity to q)
%
% To create a TWO_FLOP_JAMB_SYNC, requires
% 1. name: Jamb synchronizer name
% 2. wid: circuit width, wid(1:2)   input buffer (INV4)
%                        wid(3:4)   output buffer (INV4)
%                        wid(5:18)  Master Jamb flip-flop (JAMB_FF)
%                        wid(19:32) Slave Jamb flip-flop (JAMB_FF)
% 3. rlen: relative circuit length, use 1 by default.
% E.g. sync = TWO_FLOP_JAMB_SYNC('jlSync',450e-7,1)
classdef TWO_FLOP_JAMB_SYNC < circuit

    properties (GetAccess = 'public', SetAccess = 'private')
        vdd,gnd,d,clk,clkbar,reset; q,qbar; %Input ; Output
    end
    
    methods
        function this = TWO_FLOP_JAMB_SYNC(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1)
                widBufIn = [1;1]*wid;
                widBufOut = [1;1]*wid;
                widMaster = wid;
                widSlave = wid;
                
            else
                assert(length(wid) == 32)
                widBufIn  = wid(1:2);
                widBufOut = wid(3:4);
                widMaster = wid(5:18);
                widSlave  = wid(19:end);
                
            end
            this = this@circuit(name);
            %set the relative lengths by default to be be the same
            rlenBufferIn  = rlen;
            rlenBufferOut = rlen;
            rlenMaster = rlen;
            rlenSlave = rlen;
            
            if(iscell(rlen))
                rlenBufferIn = rlen{1};
                rlenBufferOut = rlen{2};
                rlenMaster = rlen{3};
                rlenSlave = rlen{4};
            end
            %define circuit nodes
            %inputs
            this.d = this.add_port(node('d'));
            this.clk = this.add_port(node('clk'));
            this.clkbar = this.add_port(node('clkbar'));
            this.vdd = this.add_port(node('vdd'));
            this.gnd = this.add_port(node('gnd'));
            this.reset = this.add_port(node('reset'));
            
            %internal
            dbar = this.add_port(node('dbar'));
            
            %outpus
            this.q = this.add_port(node('q'));
            this.qbar = this.add_port(node('qbar'));
            
            bufInvIn      = INV4('BufIn', widBufIn,rlenBufferIn);   this.add_element(bufInvIn);
            bufInvOut     = INV4('bufOut',widBufOut,rlenBufferOut); this.add_element(bufInvOut);
            jamb_ffMaster = JAMB_FF('JFF_0',widMaster,rlenMaster); this.add_element(jamb_ffMaster);
            jamb_ffSlave  = JAMB_FF('JFF_1',widSlave,rlenSlave); this.add_element(jamb_ffSlave);
            
            
            this.connect(this.d,bufInvIn.i);
            this.connect(dbar,bufInvIn.o,jamb_ffMaster.d);
            this.connect(jamb_ffMaster.q,jamb_ffSlave.d);
            this.connect(this.q,jamb_ffSlave.q,bufInvOut.i);
            this.connect(this.qbar,bufInvOut.o);
            this.connect(this.clk,jamb_ffMaster.clk,jamb_ffSlave.clk);
            this.connect(this.clkbar,jamb_ffMaster.clkbar,jamb_ffSlave.clkbar);
            this.connect(this.vdd,bufInvIn.vdd,bufInvOut.vdd,jamb_ffMaster.vdd,jamb_ffSlave.vdd);
            this.connect(this.gnd,bufInvIn.gnd,bufInvOut.gnd,jamb_ffMaster.gnd,jamb_ffSlave.gnd);
            this.connect(this.reset,jamb_ffMaster.reset,jamb_ffSlave.reset);
            
            this.finalize;
            
        end
        
    end
end

