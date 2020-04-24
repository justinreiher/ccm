% Definition of a STRONG_ARM_LATCH with nmos and pmos elements.
%
% The strong latch has 5 input nodes:
% 1. vdd:    power supply source
% 2. gnd:    circuit ground
% 3. d:      data input
% 4. dbar:   logical inverse of d
% 4. clk:    clock input
%
% The jamb latch has 2 output nodes:
% 1. s:      set
% 2. r:      reset (i.e. the opposite polarity to set)
%
% To create a STRONG_ARM_LATCH, requires
% 1. name: strong arm latch name
% 2. wid: circuit width, wid(1)     is the nmos device for the forward
%                                   cross-coupled inverter (nmos)
%                        wid(2)     is the pmos device for the forward
%                                   cross-coupled inverter (pmos)
%                        wid(3)     is the nmos device for the backward
%                                   cross-coupled inverter (nmos)
%                        wid(4)     is the pmos device for the backward
%                                   cross-coupled inverter (pmos)
%                        wid(5)     is the nmos device for the input signal
%                                   d nD (nmos)
%                        wid(6)     is the nmos device for the input signal
%                                   dbar nDbar (nmos)
%                        wid(7)     is the nmos pull down device for the 
%                                   clock signal nClk (nmos)
%                        wid(8)     is the nmos device whihc bridges the
%                                   the pull down path between d and dbar
%                                   nBridge (nmos)
%                        wid(9)     is the pmos device for pre-charging the
%                                   d path of the latch p0 (pmos)
%                        wid(10)    is the pmos device for pre-charging the
%                                   dbarpath of the latch p1 (pmos)
% 3. rlen: relative circuit length, use 1 by default.
% E.g. saL = STRONG_ARM_Latch('sal_0',450e-7,1) creates a strong arm latch 
% where all transistor widths are 450nm wide.
classdef STRONG_ARM_LATCH < circuit
    
    properties (GetAccess = 'public', SetAccess = 'private')
        vdd,gnd,d,dbar,clk,clkbar; s,r; %Input ; Output
    end
    
    methods
        function this = STRONG_ARM_LATCH(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1)
                
                widN_CCfwd       = wid;
                widP_CCfwd       = wid;
                widN_CCbck       = wid;
                widP_CCbck       = wid;

                %n devices
                widN_D           = wid;
                widN_Dbar        = wid;
                widN_Clk         = wid;
                widN_Bridge      = wid;
                %p devices
                widP_preCharge0 = wid;
                widP_preCharge1 = wid;
            else
                assert(length(wid) == 10)
                
                %inverter CC Fwd
                widN_CCfwd       = wid(1);
                widP_CCfwd       = wid(2);
                
                %inverter CC Bck
                widN_CCbck       = wid(3);
                widP_CCbck       = wid(4);
                
                %n devices
                widN_D           = wid(5);
                widN_Dbar        = wid(6);
                widN_Clk         = wid(7);
                widN_Bridge      = wid(8);
                %p devices
                widP_preCharge0 = wid(9);
                widP_preCharge1 = wid(10);

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
            this.s = this.add_port(node('s'));
            this.r = this.add_port(node('r'));
            
            %internal node
            z = this.add_port(node('z'));
            x = this.add_port(node('x'));
            y = this.add_port(node('y'));

            %cross-coupled inverters
            txN_CCfwd = nmos('txN_CC_Fwd','wid',widN_CCfwd,'rlen',1); this.add_element(txN_CCfwd); 
			txP_CCfwd = pmos('txP_CC_Fwd','wid',widP_CCfwd,'rlen',1); this.add_element(txP_CCfwd);
            
            txN_CCbck = nmos('txN_CC_Bck','wid',widN_CCbck,'rlen',1); this.add_element(txN_CCbck);
            txP_CCbck = pmos('txP_CC_Bck','wid',widP_CCbck,'rlen',1); this.add_element(txP_CCbck);
            
			%n devices
			txN_d 	 = nmos('txN_d','wid',widN_D,'rlen',1);	this.add_element(txN_d);
			txN_dbar	 = nmos('txN_dbar','wid',widN_Dbar,'rlen',1);	this.add_element(txN_dbar);
			txN_clk 	 = nmos('txN_clk','wid',widN_Clk,'rlen',1);	this.add_element(txN_clk);
			txN_bridge	 = nmos('txN_bridge','wid',widN_Bridge,'rlen',1);	this.add_element(txN_bridge);
			
            %p devices
            txP_preCharge0  = pmos('txP_preCharge0','wid',widP_preCharge0,'rlen',1); this.add_element(txP_preCharge0);
			txP_preCharge1	 = pmos('txP_preCharge1','wid',widP_preCharge1,'rlen',1);	this.add_element(txP_preCharge1);

            this.connect(this.s,txN_CCbck.d,txP_CCbck.d,txN_CCfwd.g,txP_CCfwd.g,txP_preCharge0.d);
            this.connect(this.r,txN_CCfwd.d,txP_CCfwd.d,txN_CCbck.g,txP_CCbck.g,txP_preCharge1.d);
            
            this.connect(x,txN_CCbck.s,txN_d.d,txN_bridge.d);
            this.connect(y,txN_CCfwd.s,txN_dbar.d,txN_bridge.s);
            this.connect(z,txN_d.s,txN_dbar.s,txN_clk.d);
            

            this.connect(this.d,txN_d.g);
			this.connect(this.dbar,txN_dbar.g);
			this.connect(this.clk,txN_clk.g,txP_preCharge0.g,txP_preCharge1.g);
			
			
            this.connect(this.vdd,txP_preCharge0.s,txP_preCharge0.b,txP_preCharge1.s,txP_preCharge1.b,txN_bridge.g,txP_CCfwd.b,txP_CCfwd.s,txP_CCbck.b,txP_CCbck.s);
            this.connect(this.gnd,txN_clk.s,txN_clk.b,txN_d.b,txN_dbar.b,txN_bridge.b,txN_CCfwd.b,txN_CCbck.b);
            
            this.finalize;
            
        end
        
    end
end

