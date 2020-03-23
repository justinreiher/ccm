classdef PASSGATE_LATCH < circuit
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
        function this = PASSGATE_LATCH(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1)
                widTot = [1,1]*wid;
                widPG0 = widTot;
                widPG1 = widTot;
                widInv0 = widTot;
                widInv1 = widTot;
                widInv2 = widTot;
            else
                assert(length(wid) == 10)
                widPG0  = wid(1:2);
                widPG1  = wid(3:4);
                widInv0 = wid(5:6);
                widInv1 = wid(7:8);
                widInv2 = wid(9:10);
            end
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
            
            %internal nodes
            x0 = this.add_port(node('x0'));
            y0 = this.add_port(node('y0'));
            
            pg_0 = PASSGATE('PG_0',widPG0,rlen); this.add_element(pg_0);
            pg_1 = PASSGATE('PG_1',widPG1,rlen); this.add_element(pg_1);
            
            inv0 = INV4('invQbar',widInv0,rlen); this.add_element(inv0);
            inv1 = INV4('fwr',widInv1,rlen); this.add_element(inv1);
            inv2 = INV4('bck',widInv2,rlen); this.add_element(inv2);

            this.connect(this.d,pg_0.i);
            this.connect(x0,pg_0.o,pg_1.o,inv1.i);
            this.connect(this.clk,pg_0.gbar,pg_1.g);
            this.connect(this.clkbar,pg_0.g,pg_1.gbar);
            this.connect(this.q,inv1.o,inv2.i,inv0.i);
            this.connect(y0,pg_1.i,inv2.o);
            this.connect(this.qbar,inv0.o);
            this.connect(this.vdd,inv0.vdd,inv1.vdd,inv2.vdd,pg_0.vdd,pg_1.vdd);
            this.connect(this.gnd,inv0.gnd,inv1.gnd,inv2.gnd,pg_0.gnd,pg_1.gnd);
            this.finalize;
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

