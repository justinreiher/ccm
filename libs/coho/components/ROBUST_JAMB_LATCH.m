classdef ROBUST_JAMB_LATCH < circuit
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
        vdd,gnd,d,clk,reset; q,qbar; %Input ; Output
    end
    
    methods
        function this = ROBUST_JAMB_LATCH(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1)
                
                widInvCCFwd  = [1;1]*wid;
                widInvCCBck  = [1;1]*wid;
                widInvOut    = [1;1]*wid;
                widNand0 = [1;1;1;1]*wid;
                widDinN     = wid;
                widClkN     = wid;
                widResetN   = wid;
                widFilterN0 = wid;
                widFilterN1 = wid;
                widBoostP0  = wid;
                widBoostP1  = wid;
                widFilterP0 = wid;
                widFilterP1 = wid;
                
            else
                assert(length(wid) == 19)
                widInvCCFwd  = wid(1:2);
                widInvCCBck  = wid(3:4);
                widInvOut    = wid(5:6);
                widNand0     = wid(7:10);
                widDinN      = wid(11);
                widClkN      = wid(12);
                widResetN    = wid(13);
                widFilterN0  = wid(14);
                widFilterN1  = wid(15);
                widBoostP0   = wid(16);
                widBoostP1   = wid(17);
                widFilterP0  = wid(18);
                widFilterP1  = wid(19);
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
            x       = this.add_port(node('x'));
            y       = this.add_port(node('y'));
            z       = this.add_port(node('z'));
            f_in    = this.add_port(node('f_in'));
            f_out   = this.add_port(node('f_out'));
            
            invCCfwd    = INV4('INV_{CC_{fwd}}',widInvCCFwd,rlen);    this.add_element(invCCfwd);
            invCCbck    = INV4('INV_{CC_{bck}}',widInvCCBck,rlen);    this.add_element(invCCbck);
            invOut      = INV4('INV_{Out}', widInvOut,rlen);          this.add_element(invOut);
            nand0       = NAND2('NAND_0',widNand0,rlen);              this.add_element(nand0);
            txN_D       = nmos('txN_{D_{in}}','wid',widDinN,'rlen',rlen);          this.add_element(txN_D);
            txN_Clk     = nmos('txN_{Clk}','wid',widClkN,'rlen',rlen);            this.add_element(txN_Clk);
            txN_Reset   = nmos('txN_{Reset}','wid',widResetN,'rlen',rlen);        this.add_element(txN_Reset);
            txN_Filter0 = nmos('txN_{Filter_0}','wid',widFilterN0,'rlen',rlen);   this.add_element(txN_Filter0);
            txN_Filter1 = nmos('txN_{Filter_1}','wid',widFilterN1,'rlen',rlen);   this.add_element(txN_Filter1);
            txP_Boost0  = pmos('txP_{Boost_0}','wid',widBoostP0,'rlen',rlen);     this.add_element(txP_Boost0);
            txP_Boost1  = pmos('txP_{Boost_1}','wid',widBoostP1,'rlen',rlen);     this.add_element(txP_Boost1);
            txP_Filter0 = pmos('txP_{Filter_0}','wid',widFilterP0,'rlen',rlen);   this.add_element(txP_Filter0);
            txP_Filter1 = pmos('txP_{Filter_1}','wid',widFilterP1,'rlen',rlen);   this.add_element(txP_Filter1);
                 
            
            this.connect(this.vdd,invCCfwd.vdd,invCCbck.vdd,invOut.vdd,nand0.vdd,...
                        txP_Boost0.b,txP_Boost1.b,txP_Filter0.b,txP_Filter1.b,...
                            txP_Boost0.s,txP_Boost1.s,txP_Filter0.s,txP_Filter1.s);
            
            this.connect(this.gnd,invCCfwd.gnd,invCCbck.gnd,invOut.gnd,nand0.gnd,...
                         txN_D.b,txN_Clk.b,txN_Reset.b,txN_Filter0.b,txN_Filter1.b,...
                         txN_Clk.s,txN_Reset.s);
            
            this.connect(this.d,txN_D.g);
            this.connect(this.clk,txN_Clk.g);
            this.connect(this.reset,txN_Reset.g);
            
            this.connect(x,invCCfwd.i,invCCbck.o,txN_D.d,txP_Boost0.d,txP_Filter0.g,txN_Filter0.g,txN_Filter1.s);
            this.connect(y,invCCfwd.o,invCCbck.i,txN_Reset.d,txP_Boost1.d,txP_Filter1.g,txN_Filter1.g,txN_Filter0.s);
            this.connect(z,txN_D.s,txN_Clk.d);
            this.connect(f_in,txN_Filter1.d,txP_Filter1.d,nand0.i2);
            this.connect(this.qbar,txN_Filter0.d,txP_Filter0.d,nand0.i1,invOut.i);
            this.connect(this.q,invOut.o);
            this.connect(f_out,nand0.o,txP_Boost0.g,txP_Boost1.g);
            
            this.finalize;
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

