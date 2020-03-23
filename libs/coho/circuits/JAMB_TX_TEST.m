classdef JAMB_TX_TEST < circuit
    % Definition to test current flowing through 2 NMOS transistors in
    % series, used to solve for the node between them
    %       vtop
    %       _|
    % d   -|_  <- txD
    %        |_ z
    %       _|
    % clk -|_  <-txClk
    %        |
    %       gnd
    % This circuit has 3 inputs:
    % 1. vdd:    power supply source, d and clk are tied to vdd
    % 2. gnd:    circuit ground
    % 3. vtop:   input from inverter pre-computed
    %
    % To create a JAMB_TX_TEST, requires
    % 1. name: circuit name
    % 2. wid: circuit width, wid(1) is for txD
    %                        wid(2) is for txClk
    % 3. rlen: relative circuit length, use 1 by default.
    
    properties (GetAccess = 'public', SetAccess = 'private')
        vdd; vtop; gnd; z;
    end
    
    methods
        function this = JAMB_TX_TEST(name,wid,rlen)
            if(nargin<1)
                error('must provide name and width');
            end
            
            if(numel(wid) == 1)
                widTxD = wid;
                widTxClk = wid;
            else
                assert(length(wid) == 2)
                widTxD = wid(1);
                widTxClk = wid(2);
            end
            
           
            this = this@circuit(name);
            
            %define circuit nodes
            %inputs
            this.vtop = this.add_port(node('vTop'));
            this.vdd  = this.add_port(node('vdd'));
            this.gnd = this.add_port(node('gnd'));
            
            %output
            this.z = this.add_port(node('z'));
          
            txD = nmos('tx1','wid',widTxD,'rlen',rlen); this.add_element(txD);
            txClk = nmos('tx2','wid',widTxClk,'rlen',rlen); this.add_element(txClk);
            
            this.connect(this.gnd,txClk.s,txD.b,txClk.b);
            this.connect(this.vtop,txD.d);
            this.connect(this.z,txD.s,txClk.d);
            this.connect(this.vdd,txClk.g,txD.g);
            
            this.finalize;
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

