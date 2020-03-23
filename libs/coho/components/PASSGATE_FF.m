classdef PASSGATE_FF < circuit
    % Definition of a PASSGATE_LATCH with PASSGATE elements and INVERTERS
    %
    % The passgate latch has 5 input nodes:
    % 1. vdd:    power supply source
    % 2. gnd:    circuit ground
    % 3. d:      data input
    % 4. clk:    clock input
    % 5. clkbar: clock bar (i.e. the opposite polarity to clock)
    % The passgate latch has 2 output nodes:
    % 1. q:      output
    % 2. qbar:   output bar (i.e. the opposite polarity to q)
    %
    % To create a PASSGATE_LATCH, requires
    % 1. name: passgate latch name
    % 2. wid: circuit width, wid(1) is for the INV
    %                        wid(2) is for the PASSGATE
    % 3. rlen: relative circuit length, use 1 by default.
    % E.g. pgLatch = PASSGATE_LATCH('pg0',[
    
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
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

