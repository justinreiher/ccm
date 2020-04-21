classdef STRONG_ARM_FF < circuit
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
        vdd,gnd,d,dbar,clk; q,qbar; %Input ; Output
    end
    
    methods
        function this = STRONG_ARM_FF(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1)
                widStrong = wid;
                widNAND0  = wid;
                widNAND1  = wid;
            else
                assert(length(wid) == 18)
                widStrong = wid(1:10);
                widNAND0  = wid(11:14);
                widNAND1  = wid(15:end);
            end
            this = this@circuit(name);
            
            %define circuit nodes
            %inputs
            this.d = this.add_port(node('d'));
            this.dbar = this.add_port(node('dbar'));
            this.clk = this.add_port(node('clk'));
            this.vdd = this.add_port(node('vdd'));
            this.gnd = this.add_port(node('gnd'));
            %outpus
            this.q = this.add_port(node('q'));
            this.qbar = this.add_port(node('qbar'));
            
            %internal
            s = this.add_port(node('s'));
            r = this.add_port(node('r'));
            
            sarmL = STRONG_ARM_LATCH('SAL',widStrong,rlen); this.add_element(sarmL);
            nand0 = NAND2('NAND_0',widNAND0,rlen); this.add_element(nand0);
            nand1 = NAND2('NAND_1',widNAND1,rlen); this.add_element(nand1);


            this.connect(this.d,sarmL.d);
            this.connect(this.dbar,sarmL.dbar);

            this.connect(s,nand0.i1,sarmL.s);
            this.connect(r,nand1.i1,sarmL.r);
            this.connect(this.clk,sarmL.clk);
            
            this.connect(this.q,nand0.o,nand1.i2);
            this.connect(this.qbar,nand1.o,nand0.i2);

            this.connect(this.vdd,sarmL.vdd,nand0.vdd,nand1.vdd);
            this.connect(this.gnd,sarmL.gnd,nand0.gnd,nand1.gnd);
            this.finalize;
            
        end
        
    end
end

