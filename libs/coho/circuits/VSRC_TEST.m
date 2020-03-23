classdef VSRC_TEST < circuit
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
        gnd; out;
    end
    
    methods
        function this = VSRC_TEST(name)
            if(nargin<1)
                error('must provide name and width');
            end
           
            this = this@circuit(name);
            
            %define circuit nodes
            %inputs
            
            this.gnd = this.add_port(node('gnd'));
            
            %outpus
            this.out = this.add_port(node('out'));
            
            
            v1 = vsrc('v1',0.1,false); this.add_element(v1);
            v2 = vsrc('v2',0.3,false); this.add_element(v2);
            v3 = vsrc('v3',0.5,false); this.add_element(v3);
            
            this.connect(this.gnd,v1.m);
            this.connect(v1.p,v2.m);
            this.connect(v2.p,v3.m);
            this.connect(this.out,v3.p);
            
            this.finalize;
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

