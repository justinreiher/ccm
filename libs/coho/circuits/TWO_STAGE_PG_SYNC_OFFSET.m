classdef TWO_STAGE_PG_SYNC_OFFSET < circuit
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
        function this = TWO_STAGE_PG_SYNC_OFFSET(name,offset,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1), wid=[1;1]*wid; end
            this = this@circuit(name);
            
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
            
            
            pg_ffMaster = PASSGATE_FF_OFFSET('PGFF_0',offset,wid(2),rlen); this.add_element(pg_ffMaster);
            pg_ffSlave  = PASSGATE_FF('PGFF_1',wid(2),rlen); this.add_element(pg_ffSlave);
            
            this.connect(this.d,pg_ffMaster.d);
            this.connect(pg_ffMaster.q,pg_ffSlave.d);
            this.connect(this.clk,pg_ffMaster.clk,pg_ffSlave.clk);
            this.connect(this.clkbar,pg_ffMaster.clkbar,pg_ffSlave.clkbar);
            this.connect(this.vdd,pg_ffMaster.vdd,pg_ffSlave.vdd);
            this.connect(this.gnd,pg_ffMaster.gnd,pg_ffSlave.gnd);
            
            
            this.connect(this.q,pg_ffSlave.q);
            this.connect(this.qbar,pg_ffSlave.qbar);
            this.finalize;
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

