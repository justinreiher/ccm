classdef TX_TEST_TWO < circuit
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
        vdd; gnd; in; in1; in2; out;
    end
    
    methods
        function this = TX_TEST_TWO(name,wid,rlen)
            if(nargin<1)
                error('must provide name and width');
            end
            
            if(numel(wid) == 1)
                widInv0 = [1,1]*wid;
                widInv1 = [1,1]*wid;
                widTx1 = wid;
                widTx2 = wid;
            else
                assert(length(wid) == 4)
                widInv0 = wid(1:2);
                widInv1 = wid(3:4);
                widTx1 = wid(3);
                widTx2 = wid(4);
            end
            
           
            this = this@circuit(name);
            
            %define circuit nodes
            %inputs
            this.vdd = this.add_port(node('vdd'));
            this.gnd = this.add_port(node('gnd'));
            this.in  = this.add_port(node('in'));
            this.in1 = this.add_port(node('in1'));
            this.in2 = this.add_port(node('in2'));
            
            %internal nodes
            x0 = this.add_port(node('x0'));
            y0 = this.add_port(node('y0'));
            %outpus
            this.out = this.add_port(node('out'));
            
            
            inv0 = INV4('inv0',widInv0,rlen); this.add_element(inv0);
            inv1 = INV4('inv1',widInv1,rlen); this.add_element(inv1);
            tx1 = nmos('tx1','wid',widTx1,'rlen',rlen); this.add_element(tx1);
            tx2 = nmos('tx2','wid',widTx2,'rlen',rlen); this.add_element(tx2);
            
            this.connect(this.gnd,tx2.s,tx1.b,tx2.b,inv0.gnd,inv1.gnd);
            this.connect(this.vdd,inv0.vdd,inv1.vdd);
            this.connect(this.in1,tx1.g);
            this.connect(this.in2,tx2.g);
            this.connect(this.in,inv0.i);
            this.connect(x0,tx1.s,tx2.d,inv1.i);
            this.connect(y0,inv0.o,tx1.d);
            
            this.connect(this.out,inv1.o);
            
            this.finalize;
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

