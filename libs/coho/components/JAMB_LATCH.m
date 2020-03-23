classdef JAMB_LATCH < circuit
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
        vdd,gnd,d,reset,clk; q,qbar; %Input ; Output
    end
    
    methods
        function this = JAMB_LATCH(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1)
                widTot = [1,1]*wid;
                widInvCCFwr  = widTot;
                widInvCCBck  = widTot;
                wtxDin   = wid;
                wtxClk   = wid;
                wtxReset = wid;
                
            else
                assert(length(wid) == 7)
                widInvCCFwr  = wid(1:2);
                widInvCCBck  = wid(3:4);
                wtxDin   = wid(5);
                wtxClk   = wid(6);
                wtxReset = wid(7);
            end
            this = this@circuit(name);
            
            %define circuit nodes
            %inputs
            this.d = this.add_port(node('d'));
            this.clk = this.add_port(node('clk'));
            this.vdd = this.add_port(node('vdd'));
            this.gnd = this.add_port(node('gnd'));
            this.reset = this.add_port(node('reset'));
            
            %outpus
            this.q = this.add_port(node('q'));
            this.qbar = this.add_port(node('qbar'));
            
            %internal node
            z0   = this.add_port(node('z0'));
            
            inv0 = INV4('Fwr',widInvCCFwr,rlen); this.add_element(inv0);
            inv1 = INV4('Bck',widInvCCBck,rlen); this.add_element(inv1);
            txDin  = nmos('tx_D','wid',wtxDin,'rlen',1); this.add_element(txDin);
            txClk = nmos('tx_Clk','wid',wtxClk,'rlen',1); this.add_element(txClk);
            txReset  = nmos('tx_Reset','wid',wtxReset,'rlen',1); this.add_element(txReset);

            this.connect(this.d,txDin.g);
            this.connect(this.reset,txReset.g);
            this.connect(this.qbar,inv0.i,inv1.o,txDin.d);
            this.connect(this.q,inv0.o,inv1.i,txReset.d);
            this.connect(this.clk,txClk.g);
            this.connect(z0,txDin.s,txClk.d);
           
            this.connect(this.vdd,inv0.vdd,inv1.vdd);
            this.connect(this.gnd,inv0.gnd,inv1.gnd,txDin.b,txReset.b,txClk.b,txClk.s,txReset.s);
            
            this.finalize;
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

