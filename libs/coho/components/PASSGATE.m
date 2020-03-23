classdef PASSGATE < circuit
    %Circuit description for a Passgate
    
    properties (GetAccess='public', SetAccess='private')
        vdd,gnd,i,g,gbar; o; %Input ; Output
    end
    
    methods
        function this = PASSGATE(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
             error('must provide name and width'); 
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1), wid = [1;1]*wid; end
            this = this@circuit(name); %super-class constructor
           
            %define circuit nodes            
            %inputs
            this.i = this.add_port(node('i'));
            this.g = this.add_port(node('g'));
            this.gbar = this.add_port(node('gbar'));
            this.vdd = this.add_port(node('vdd'));
            this.gnd = this.add_port(node('gnd'));
            %outputs
            this.o = this.add_port(node('q'));
            
            %create devices wid is of the form wid = [widN,widP] if wid is
            %called as a single number then both the P and N device will
            %have the same width.
            
            widN = wid(1);
            widP = wid(2);
            
            n1 = nmos('n1','wid',widN,'rlen',rlen); this.add_element(n1);
            p1 = pmos('p1','wid',widP,'rlen',rlen); this.add_element(p1);
            
            %define connections
            this.connect(this.i,n1.s,p1.s);
            this.connect(this.g,n1.g);
            this.connect(this.gbar,p1.g);
            this.connect(this.o,n1.d,p1.d);
            this.connect(this.vdd,p1.b);
            this.connect(this.gnd,n1.b);
            
            this.finalize; % complete component
            
        end
    end
end

