% Definition of a TWO_FLOP_ROBUST_JAMB_SYNC with ROBUST_JAMB_FF and INV4 elements
% and defines a two flip-flop Robust Jamb latch synchronizers 
% The Robust Jamb synchronizer has 6 input nodes:
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
% To create a TWO_FLOP_ROBUST_JAMB_SYNC, requires
% 1. name: Jamb synchronizer name
% 2. wid: circuit width, wid(1:2)   input buffer (INV4)
%                        wid(3:40)  Master Jamb flip-flop (ROBUST_JAMB_FF)
%                        wid(41:78) Slave Jamb flip-flop (ROBUST_JAMB_FF)
% 3. rlen: relative circuit length, use 1 by default.
% E.g. sync = TWO_FLOP_ROBUST_JAMB_SYNC('rjlSync',450e-7,1)
classdef TWO_FLOP_ROBUST_JAMB_SYNC < circuit
    
    properties (GetAccess = 'public', SetAccess = 'private')
        vdd,gnd,d,clk,clkbar,reset; q,qbar; %Input ; Output
    end
    
    methods
        function this = TWO_FLOP_ROBUST_JAMB_SYNC(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1)
                widBuf = [1;1]*wid;
                widMaster = wid;
                widSlave = wid;
            else
                assert(length(wid) == 78)
                widBuf = wid(1:2);
                widMaster = wid(3:40);
                widSlave  = wid(41:end);
            end
            this = this@circuit(name);
            %set the relative lengths by default to be be the same
            rlenMaster = rlen;
            rlenSlave = rlen;
            rlenBuf = rlen;
            if(iscell(rlen))
                rlenMaster = rlen{1};
                rlenSlave = rlen{2};
                rlenBuf  = rlen{3};
            end
            %define circuit nodes
            %inputs
            this.d = this.add_port(node('d'));
            this.clk = this.add_port(node('clk'));
            this.clkbar = this.add_port(node('clkbar'));
            this.reset = this.add_port(node('reset'));
            this.vdd = this.add_port(node('vdd'));
            this.gnd = this.add_port(node('gnd'));
            
            %outpus
            this.q = this.add_port(node('q'));
            this.qbar = this.add_port(node('qbar'));
            
            %internal
            dbar = this.add_port(node('dbar'));
            
            invBuf        = INV4('inv_{Buf}',widBuf,rlenBuf); this.add_element(invBuf);
            jamb_ffMaster = ROBUST_JAMB_FF('RJFF_0',widMaster,rlenMaster); this.add_element(jamb_ffMaster);
            jamb_ffSlave  = ROBUST_JAMB_FF('RJFF_1',widSlave,rlenSlave); this.add_element(jamb_ffSlave);
            
            this.connect(this.d,invBuf.i);
            this.connect(dbar,invBuf.o,jamb_ffMaster.d);
            this.connect(jamb_ffMaster.qbar,jamb_ffSlave.d);

            this.connect(this.clk,jamb_ffMaster.clk,jamb_ffSlave.clk);
            this.connect(this.clkbar,jamb_ffMaster.clkbar,jamb_ffSlave.clkbar);
            this.connect(this.reset,jamb_ffMaster.reset,jamb_ffSlave.reset);
            
            this.connect(this.vdd,jamb_ffMaster.vdd,jamb_ffSlave.vdd,invBuf.vdd);
            this.connect(this.gnd,jamb_ffMaster.gnd,jamb_ffSlave.gnd,invBuf.gnd);
            
            
            this.connect(this.qbar,jamb_ffSlave.qbar);
            this.connect(this.q,jamb_ffSlave.q);
            this.finalize;
            
        end
        
    end
end

