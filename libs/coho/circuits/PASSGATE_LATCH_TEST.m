classdef PASSGATE_LATCH_TEST < circuit
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
        function this = PASSGATE_LATCH_TEST(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1)
                widBuf = [1,1]*wid;
                widLatch = wid;
            else
                assert(length(wid) == 12)
                widBuf  = wid(1:2);
                widLatch = wid(3:end);
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
            
            %internal node
            dbar = this.add_port(node('dbar'));
           
            
            invBuf = INV4('buf',widBuf,rlen); this.add_element(invBuf);
            pgl    = PASSGATE_LATCH('pgl0',widLatch,rlen); this.add_element(pgl);
            

            this.connect(this.d,invBuf.i);
            this.connect(dbar,invBuf.o,pgl.d);
            this.connect(this.clk,pgl.clk);
            this.connect(this.clkbar,pgl.clkbar);
            this.connect(this.vdd,invBuf.vdd,pgl.vdd);
            this.connect(this.gnd,invBuf.gnd,pgl.gnd);
            
            this.connect(this.q,pgl.q);
            this.connect(this.qbar,pgl.qbar);
            this.finalize;
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

