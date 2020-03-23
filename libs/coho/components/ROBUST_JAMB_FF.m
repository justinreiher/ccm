classdef ROBUST_JAMB_FF < circuit
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
        vdd,gnd,d,clk,clkbar,reset; q,qbar; %Input ; Output
    end
    
    methods
        function this = ROBUST_JAMB_FF(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1)
                widMaster = wid;
                widSlave  = wid;            
            else
                assert(length(wid) == 38)
                widMaster = wid(1:19);
                widSlave  = wid(20:end);
            end
            this = this@circuit(name);
            
            %define circuit nodes
            %inputs
            this.d = this.add_port(node('d'));
            this.clk = this.add_port(node('clk'));
            this.clkbar = this.add_port(node('clkbar'));
            this.reset  = this.add_port(node('reset'));
            this.vdd = this.add_port(node('vdd'));
            this.gnd = this.add_port(node('gnd'));
            
            %outpus
            this.q = this.add_port(node('q'));
            this.qbar = this.add_port(node('qbar'));
            
            jl_master = ROBUST_JAMB_LATCH('RJL_0',widMaster,rlen); this.add_element(jl_master);
            jl_slave = ROBUST_JAMB_LATCH('RJL_1',widSlave,rlen); this.add_element(jl_slave);

            this.connect(this.d,jl_master.d);
            this.connect(jl_master.qbar,jl_slave.d);
            this.connect(this.qbar,jl_slave.qbar);
            this.connect(this.q,jl_slave.q);
            this.connect(this.clk,jl_master.clk);
            this.connect(this.clkbar,jl_slave.clk);
            this.connect(this.reset,jl_slave.reset,jl_master.reset);

            this.connect(this.vdd,jl_master.vdd,jl_slave.vdd);
            this.connect(this.gnd,jl_master.gnd,jl_slave.gnd);
            this.finalize;
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

