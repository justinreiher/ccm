% Definition of a JAMB_FF_BUFFERED with JAMB_LATCH and INV4 elements
% defines a buffered Jamb flip-flop where there are inverters at the input
% to each stage to buffer them.
%
% The Jamb flip-flop has 6 input nodes:
% 1. vdd:    power supply source
% 2. gnd:    circuit ground
% 3. d:      data input
% 4. clk:    clock input
% 5. clkbar: clock bar (i.e. the opposite polarity to clock)
% 6. reset:  the reset signal to reset the Jamb latches
% The passgate latch has 2 output nodes:
% 1. q:      output
% 2. qbar:   output bar (i.e. the opposite polarity to q)
%
% To create a JAMB_FF_BUFFERED, requires
% 1. name: Jamb flip-flop name
% 2. wid: circuit width, wid(1:2)   is for the input buffer (INV4)
%                        wid(3:4)   is for the buffer between Master and
%                                   Slave (INV4)
%                        wid(5:11)  is for the Master stage (JAMB_LATCH)
%                        wid(12:18) is for the Slave stage (JAMB_LATCH)
% 3. rlen: relative circuit length, use 1 by default.
% E.g. bufJLFF = JAMB_FF_BUFFERED('bufJLFF_0',450e-7,1);
classdef JAMB_FF_BUFFERED < circuit
    
    properties (GetAccess = 'public', SetAccess = 'private')
        vdd,gnd,d,clk,clkbar,reset; q,qbar; %Input ; Output
    end
    
    methods
        function this = JAMB_FF_BUFFERED(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1)
                widMaster = wid;
                widSlave  = wid;
                
            else
                assert(length(wid) == 18)
                widBufIn  = wid(1:2);
                widBufMid = wid(3:4);
                widMaster = wid(5:11);
                widSlave  = wid(12:end);
            end
            this = this@circuit(name);
            
            %define circuit nodes
            %inputs
            this.d = this.add_port(node('d'));
            this.clk = this.add_port(node('clk'));
            this.clkbar = this.add_port(node('clkbar'));
            this.vdd = this.add_port(node('vdd'));
            this.gnd = this.add_port(node('gnd'));
            this.reset = this.add_port(node('reset'));
            dbar = this.add_port(node('dbar'));
            
            %outpus
            this.q = this.add_port(node('q'));
            this.qbar = this.add_port(node('qbar'));
           
            invBufIn  = INV4('bufIn',widBufIn,rlen); this.add_element(invBufIn);
            invBufMid = INV4('bufMid',widBufMid,rlen);this.add_element(invBufMid);
            jl_master = JAMB_LATCH('JL_0',widMaster,rlen); this.add_element(jl_master);
            jl_slave = JAMB_LATCH('JL_1',widSlave,rlen); this.add_element(jl_slave);
            
            this.connect(this.d,invBufIn.i);
            this.connect(dbar,invBufIn.o,jl_master.d);
            this.connect(jl_master.qbar,invBufMid.i)
            this.connect(invBufMid.o,jl_slave.d);
            this.connect(this.qbar,jl_slave.qbar);
            this.connect(this.q,jl_slave.q);
            this.connect(this.clk,jl_master.clk);
            this.connect(this.clkbar,jl_slave.clk);

            this.connect(this.vdd,jl_master.vdd,jl_slave.vdd,invBufIn.vdd,invBufMid.vdd);
            this.connect(this.gnd,jl_master.gnd,jl_slave.gnd,invBufIn.gnd,invBufMid.gnd);
            this.connect(this.reset,jl_master.reset,jl_slave.reset);
            this.finalize;
            
        end
    end
end

