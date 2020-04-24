% Definition of a PASSGATE_FF with PASSGATE_LATCH elements
%
% The passgate flip-flop has 5 input nodes:
% 1. vdd:    power supply source
% 2. gnd:    circuit ground
% 3. d:      data input
% 4. clk:    clock input
% 5. clkbar: clock bar (i.e. the opposite polarity to clock)
% The passgate latch has 2 output nodes:
% 1. q:      output
% 2. qbar:   output bar (i.e. the opposite polarity to q)
%
% To create a PASSGATE_FF, requires
% 1. name: passgate flip-flop name
% 2. wid: circuit width, wid(1:10)  is for the master stage
%                        wid(11:20) is for the slave stage
% 3. rlen: relative circuit length, use 1 by default.
% E.g. pgFF = PASSGATE_FF('pgFF_0',450e-7,1) creates a passgate flip-flop
% where all the transistor widths are 450nm.
classdef PASSGATE_FF < circuit
    
    properties (GetAccess = 'public', SetAccess = 'private')
        vdd,gnd,d,clk,clkbar; q,qbar; %Input ; Output
    end
    
    methods
        function this = PASSGATE_FF(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1)
                widMaster = wid;
                widSlave = wid;
            else
                assert(length(wid) == 20)
                widMaster = wid(1:10);
                widSlave  =wid(11:end);
            end
            this = this@circuit(name);
            
            rlenMaster = rlen;
            rlenSlave = rlen;
            
            if(iscell(rlen))
                rlenMaster = rlen{1};
                rlenSlave = rlen{2};
            end
            
            %define circuit nodes
            %inputs
            this.d = this.add_port(node('d'));
            this.clk = this.add_port(node('clk'));
            this.clkbar = this.add_port(node('clkbar'));
            this.vdd = this.add_port(node('vdd'));
            this.gnd = this.add_port(node('gnd'));
            
            %outpus
            this.q = this.add_port(node('q'));
            this.qbar = this.add_port(node('qbar'));
            
            
            pgl_master = PASSGATE_LATCH('PGL_0',widMaster,rlenMaster); this.add_element(pgl_master);
            pgl_slave = PASSGATE_LATCH('PGL_1',widSlave,rlenSlave); this.add_element(pgl_slave);
            
            this.connect(this.d,pgl_master.d);
            this.connect(pgl_master.qbar,pgl_slave.d);
            this.connect(this.clk,pgl_master.clk,pgl_slave.clkbar);
            this.connect(this.clkbar,pgl_master.clkbar,pgl_slave.clk);
            this.connect(this.vdd,pgl_master.vdd,pgl_slave.vdd);
            this.connect(this.gnd,pgl_master.gnd,pgl_slave.gnd);
            
            
            this.connect(this.q,pgl_slave.q);
            this.connect(this.qbar,pgl_slave.qbar);
            this.finalize;
            
        end
    end
end

