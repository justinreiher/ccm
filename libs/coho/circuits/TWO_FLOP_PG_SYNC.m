classdef TWO_FLOP_PG_SYNC < circuit
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
        function this = TWO_FLOP_PG_SYNC(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1)
                widBuf=[1;1]*wid;
                widMaster = wid;
                widSlave = wid;
            else
                assert(length(wid) == 42)
                widBuf = wid(1:2);
                widMaster = wid(3:22);
                widSlave  = wid(23:end);
            end
            this = this@circuit(name);
            %set the relative lengths by default to be be the same
            rlenMaster = rlen;
            rlenSlave = rlen;
            rlenBuffer = rlen;
            if(iscell(rlen))
                rlenBuffer = rlen{1};
                rlenMaster = rlen{2};
                rlenSlave = rlen{3};
            end
            %define circuit nodes
            %inputs
            this.d = this.add_port(node('d'));
            this.clk = this.add_port(node('clk'));
            this.clkbar = this.add_port(node('clkbar'));
            this.vdd = this.add_port(node('vdd'));
            this.gnd = this.add_port(node('gnd'));
            
            %internal node
            dbar = this.add_port(node('dbar'));
            
            %outpus
            this.q = this.add_port(node('q'));
            this.qbar = this.add_port(node('qbar'));
            
            bufInv = INV4('buffer',widBuf,rlenBuffer); this.add_element(bufInv);
            pg_ffMaster = PASSGATE_FF('PGFF_0',widMaster,rlenMaster); this.add_element(pg_ffMaster);
            pg_ffSlave  = PASSGATE_FF('PGFF_1',widSlave,rlenSlave); this.add_element(pg_ffSlave);
            
            this.connect(this.d,bufInv.i);
            this.connect(dbar,bufInv.o,pg_ffMaster.d);
            this.connect(pg_ffMaster.qbar,pg_ffSlave.d);
            this.connect(this.clk,pg_ffMaster.clk,pg_ffSlave.clk);
            this.connect(this.clkbar,pg_ffMaster.clkbar,pg_ffSlave.clkbar);
            this.connect(this.vdd,pg_ffMaster.vdd,pg_ffSlave.vdd,bufInv.vdd);
            this.connect(this.gnd,pg_ffMaster.gnd,pg_ffSlave.gnd,bufInv.gnd);
            
            
            this.connect(this.q,pg_ffSlave.q);
            this.connect(this.qbar,pg_ffSlave.qbar);
            this.finalize;
            
        end
        
    end
end

