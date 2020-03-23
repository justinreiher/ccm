classdef PSEUDO_NMOS_FF < circuit
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
        vref,vdd,gnd,d,dbar,clk,clkbar; q,qbar; %Input ; Output
    end
    
    methods
        function this = PSEUDO_NMOS_FF(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1)
                widMaster = wid;
                widSlave  = wid;
            else
                assert(length(wid) == 18)
                widMaster = wid(1:9);
                widSlave  = wid(10:end);
            end
            this = this@circuit(name);
            
            %define circuit nodes
            %inputs
            this.d = this.add_port(node('d'));
            this.dbar = this.add_port(node('dbar'));
            this.clk = this.add_port(node('clk'));
            this.clkbar = this.add_port(node('clkbar'));
            this.vdd = this.add_port(node('vdd'));
            this.gnd = this.add_port(node('gnd'));
            this.vref = this.add_port(node('vref'));
            %outpus
            this.q = this.add_port(node('q'));
            this.qbar = this.add_port(node('qbar'));
            
            
            pNmos_master = PSEUDO_NMOS_LATCH('pnL_M',widMaster,rlen); this.add_element(pNmos_master);
            pNmos_slave = PSEUDO_NMOS_LATCH('pnL_S',widSlave,rlen); this.add_element(pNmos_slave);


            this.connect(this.d,pNmos_master.d);
            this.connect(this.dbar,pNmos_master.dbar);
            this.connect(pNmos_master.q,pNmos_slave.d);
            this.connect(pNmos_master.qbar,pNmos_slave.dbar);
            this.connect(this.qbar,pNmos_slave.qbar);
            this.connect(this.q,pNmos_slave.q);
            this.connect(this.clkbar,pNmos_master.clk);
            this.connect(this.clk,pNmos_slave.clk);

            this.connect(this.vdd,pNmos_master.vdd,pNmos_slave.vdd);
            this.connect(this.gnd,pNmos_master.gnd,pNmos_slave.gnd);
            this.connect(this.vref,pNmos_master.vref,pNmos_slave.vref);
            this.finalize;
            
        end
        
    end
end

