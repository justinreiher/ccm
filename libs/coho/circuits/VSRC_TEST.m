% Definition of a VSRC_TEST with vrsc elements to test having voltage
% source in the middle of a circuit, making sure that the voltages add
% properly
%
% The test has 1 input node:
% 1. gnd:    circuit ground
% The test has 1 output nodes:
% 1. out:      output
%
% To create a VSRC_TEST, requires
% 1. name: test name
% E.g. myTest = VSRC_TEST('test')
classdef VSRC_TEST < circuit
    
    properties (GetAccess = 'public', SetAccess = 'private')
        gnd; out;
    end
    
    methods
        function this = VSRC_TEST(name)
            if(nargin<1)
                error('must provide name');
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
    end
end

