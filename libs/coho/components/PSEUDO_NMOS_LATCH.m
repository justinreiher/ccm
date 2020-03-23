classdef PSEUDO_NMOS_LATCH < circuit
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
        vref,vdd,gnd,d,dbar,clk; q,qbar; %Input ; Output
    end
    
    methods
        function this = PSEUDO_NMOS_LATCH(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1)
                widPull_0  = wid;
                widPull_1  = wid;
                widQ       = wid;
                widQbar    = wid;
                widD       = wid;
                widDbar    = wid;
                widClkD    = wid;
                widClkDbar = wid;
                widBridge  = wid;    
            else
                assert(length(wid) == 9)
                widPull_0  = wid(1);
                widPull_1  = wid(2);
                widQ       = wid(3);
                widQbar    = wid(4);
                widD       = wid(5);
                widDbar    = wid(6);
                widClkD    = wid(7);
                widClkDbar = wid(8);
                widBridge  = wid(9);
            end
            this = this@circuit(name);
            
            %define circuit nodes
            %inputs
            this.d = this.add_port(node('d'));
            this.dbar = this.add_port(node('dbar'));
            this.clk = this.add_port(node('clk'));
            this.vdd = this.add_port(node('vdd'));
            this.gnd = this.add_port(node('gnd'));
            this.vref = this.add_port(node('vref'));
            %outpus
            this.q = this.add_port(node('q'));
            this.qbar = this.add_port(node('qbar'));
            
            %internal node
            mm = this.add_port(node('mm'));
            mm_bar   = this.add_port(node('mm_bar'));
            
            pullUp_0  = pmos('pullup0','wid',widPull_0,'rlen',1); this.add_element(pullUp_0);
            pullUp_1  = pmos('pullup1','wid',widPull_1,'rlen',1); this.add_element(pullUp_1);
            txQ       = nmos('txQ','wid',widQ,'rlen',1); this.add_element(txQ);
            txQbar    = nmos('txQbar','wid',widQbar,'rlen',1); this.add_element(txQbar);
            txD       = nmos('txD','wid',widD,'rlen',1); this.add_element(txD);
            txDbar    = nmos('txDbar','wid',widDbar,'rlen',1); this.add_element(txDbar);
            txClkD    = nmos('txClkD','wid',widClkD,'rlen',1); this.add_element(txClkD);
            txClkDbar = nmos('txClkDbar','wid',widClkDbar,'rlen',1); this.add_element(txClkDbar);
            txBridge  = nmos('txBridge','wid',widBridge,'rlen',1); this.add_element(txBridge);
            
            %input and output connections
            this.connect(this.d,txD.g);
            this.connect(this.dbar,txDbar.g);
            this.connect(this.clk,txClkD.g,txClkDbar.g,txBridge.g);
            this.connect(this.q,txQ.g,pullUp_1.d,txQbar.d);
            this.connect(this.qbar,txQbar.g,pullUp_0.d,txQ.d);
            
            %power connections
            this.connect(this.gnd,txD.s,txDbar.s,txClkD.s,txClkDbar.s,...
                         txQ.b,txQbar.b,txD.b,txDbar.b,txDbar.b,txClkD.b,txClkDbar.b,txBridge.b);
            this.connect(this.vref,pullUp_0.g,pullUp_1.g);
            this.connect(this.vdd,pullUp_0.s,pullUp_1.s,pullUp_0.b,pullUp_1.b);
            
            %internal connections
            this.connect(mm,txBridge.d,txClkDbar.d,txQbar.s,txDbar.d);
            this.connect(mm_bar,txBridge.s,txQbar.s,txClkD.d,txD.d);
            
            this.finalize;
            
        end
        
    end
end

