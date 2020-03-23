classdef TWO_FLOP_PSEUDO_NMOS_SYNC_TEST < circuit
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
        vdd,gnd,d,dbar,clk,clkbar; q,qbar; %Input ; Output
    end
    
    methods
        function this = TWO_FLOP_PSEUDO_NMOS_SYNC_TEST(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1)
                widBufIn0=[1;1]*wid;
                widBufIn1 = [1;1]*wid;
                widMaster = wid;
                widSlave = wid;
                widBufOut0=[1;1]*wid;
                widBufOut1 = [1;1]*wid;
            else
                assert(length(wid) == 44)
                widBufIn0 = wid(1:2);
                widBufIn1 = wid(3:4);
                widMaster = wid(5:22);
                widSlave  = wid(23:40);
                widBufOut0 = wid(41:42);
                widBufOut1 = wid(42:44);

            end
            this = this@circuit(name);
            %set the relative lengths by default to be be the same
            rlenBufIn0 = rlen;
            rlenBufIn1 = rlen;
            rlenMaster = rlen;
            rlenSlave = rlen;
            rlenBufOut0 = rlen;
            rlenBufOut1 = rlen;
            if(iscell(rlen))
                rlenBufIn0 = rlen{1};
                rlenBufIn1 = rlen{2};
                rlenMaster = rlen{3};
                rlenSlave = rlen{4};
                rlenBufOut0 = rlen{5};
                rlenBufOut1 = rlen{6};
            end
            %define circuit nodes
            %inputs
            this.d = this.add_port(node('d'));
            this.dbar = this.add_port(node('dbar'));
            this.clk = this.add_port(node('clk'));
            this.clkbar = this.add_port(node('clkbar'));
            this.vdd = this.add_port(node('vdd'));
            this.gnd = this.add_port(node('gnd'));
            
            %internal node
            dinBar_buf = this.add_port(node('dinBar_buf'));
            din_buf    = this.add_port(node('din_buf'));
            %outpus
            this.q = this.add_port(node('q'));
            this.qbar = this.add_port(node('qbar'));
            
            dinBuf        = INV4('dinBuf',widBufIn0,rlenBufIn0); this.add_element(dinBuf);
            dinBarBuf     = INV4('dinBarBuf',widBufIn1,rlenBufIn1); this.add_element(dinBarBuf);
            pNmos_ffMaster = PSEUDO_NMOS_FF('PSN_0',widMaster,rlenMaster); this.add_element(pNmos_ffMaster);
            pNmos_ffSlave  = PSEUDO_NMOS_FF('PSN_1',widSlave,rlenSlave); this.add_element(pNmos_ffSlave);
            bufQOut     = INV4('bufQOut',widBufOut0,rlenBufOut0); this.add_element(bufQOut);
            bufQBarOut  = INV4('bufQBarOut',widBufOut1,rlenBufOut1); this.add_element(bufQBarOut);
            
            this.connect(this.d,dinBarBuf.i);
            this.connect(this.dbar,dinBuf.i);
            this.connect(din_buf,dinBuf.o,pNmos_ffMaster.d);
            this.connect(dinBar_buf,dinBarBuf.o,pNmos_ffMaster.dbar);
            this.connect(pNmos_ffMaster.q,pNmos_ffSlave.d);
            this.connect(pNmos_ffMaster.qbar,pNmos_ffSlave.dbar);
            this.connect(this.clk,pNmos_ffMaster.clk,pNmos_ffSlave.clk);
            this.connect(this.clkbar,pNmos_ffMaster.clkbar,pNmos_ffSlave.clkbar);
            this.connect(this.vdd,pNmos_ffMaster.vdd,pNmos_ffSlave.vdd,dinBuf.vdd,dinBarBuf.vdd,bufQout.vdd,bufQBarOut.vdd);
            this.connect(this.gnd,pNmos_ffMaster.gnd,pNmos_ffSlave.gnd,dinBuf.gnd,dinBarBuf.gnd,bufQout.gnd,bufQBarOut.gnd);
            
            this.connect(bufQout.i,pNmos_ffSlave.qbar);
            this.connect(bufQBarOut.i,pNmos_ffSlave.q);
            this.connect(this.q,bufQout.o);
            this.connect(this.qbar,bufQBarOut.o);
            this.finalize;
            
        end
        
    end
end

