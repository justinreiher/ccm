% Definition of a THREE_FLOP_PG_SYNC with PASSGATE_FF and INV4 elements
% defines a three flip-flop passgate synchronizer
%
% The three stage passgate synchronizer has 5 input nodes:
% 1. vdd:    power supply source
% 2. gnd:    circuit ground
% 3. d:      data input
% 4. clk:    clock input
% 5. clkbar: clock bar (i.e. the opposite polarity to clock)
% The synchronizer has 2 output nodes:
% 1. q:      output
% 2. qbar:   output bar (i.e. the opposite polarity to q)
%
% To create a THREE_FLOP_PG_SYNC, requires
% 1. name: synchronizer latch name
% 2. wid: circuit width, wid(1:2)   is for the buffer INV4
%                        wid(3:22)  is for the first Master PASSGATE_FF
%                        wid(23:42) is for the first Slave PASSGATE_FF
%                        wid(43:62) is for the second Master PASSGATE_FF
% 3. rlen: relative circuit length, use 1 by default.
% E.g. sync = THREE_FLOP_PG_SYNC('syncPG3',450e-7,1) creates a 3 stage
% passgate synchronizer with all transistor widths set to 450nm.

classdef THREE_FLOP_PG_SYNC < circuit
    
    properties (GetAccess = 'public', SetAccess = 'private')
        vdd,gnd,d,clk,clkbar; q,qbar; %Input ; Output
    end
    
    methods
        function this = THREE_FLOP_PG_SYNC(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1)
                widBuf=[1;1]*wid;
                widMaster1 = wid;
                widSlave1 = wid;
                widMaster2 = wid;
            else
                assert(length(wid) == 62)
                widBuf = wid(1:2);
                widMaster1 = wid(3:22);
                widSlave1  = wid(23:42);
                widMaster2 = wid(43:end);
            end
            this = this@circuit(name);
            %set the relative lengths by default to be be the same
            rlenMaster1 = rlen;
            rlenSlave1 = rlen;
            rlenMaster2 = rlen;
            rlenBuffer = rlen;
            if(iscell(rlen))
                rlenBuffer = rlen{1};
                rlenMaster1 = rlen{2};
                rlenSlave1 = rlen{3};
                rlenMaster2 = rlen{4};
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
            pg_ffMaster1 = PASSGATE_FF('PGFF_0',widMaster1,rlenMaster1); this.add_element(pg_ffMaster1);
            pg_ffSlave1  = PASSGATE_FF('PGFF_1',widSlave1,rlenSlave1); this.add_element(pg_ffSlave1);
            pg_ffMaster2 = PASSGATE_FF('PGFF_2',widMaster2,rlenMaster2); this.add_element(pg_ffMaster2);
            
            this.connect(this.d,bufInv.i);
            this.connect(dbar,bufInv.o,pg_ffMaster1.d);
            this.connect(pg_ffMaster1.qbar,pg_ffSlave1.d);
            this.connect(pg_ffMaster2.d,pg_ffSlave1.qbar);
            this.connect(this.clk,pg_ffMaster1.clk,pg_ffSlave1.clk,pg_ffMaster2.clk);
            this.connect(this.clkbar,pg_ffMaster1.clkbar,pg_ffSlave1.clkbar,pg_ffMaster2.clkbar);
            this.connect(this.vdd,pg_ffMaster1.vdd,pg_ffSlave1.vdd,pg_ffMaster2.vdd,bufInv.vdd);
            this.connect(this.gnd,pg_ffMaster1.gnd,pg_ffSlave1.gnd,pg_ffMaster2.gnd,bufInv.gnd);
            
            
            this.connect(this.q,pg_ffMaster2.q);
            this.connect(this.qbar,pg_ffMaster2.qbar);
            this.finalize;
            
        end
        
    end
end

