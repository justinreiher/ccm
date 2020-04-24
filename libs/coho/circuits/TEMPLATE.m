% Template to be used for your own circuit
classdef TEMPLATE < circuit
    
    properties (GetAccess = 'public', SetAccess = 'private')
        %Input; %outputs
    end
    
    methods
        function this = TEMPLATE(name,wid,rlen)
            if(nargin<2||isempty(name)||isempty(wid))
                error('must provide name and width');
            end
            if(nargin<3||isempty(rlen)), rlen = 1; end
            if(numel(wid)==1)
            else
                assert(length(wid) == 0) %replace 0 with the correct number
            end
            this = this@circuit(name);
            
            this.finalize;
            
        end
        
    end
end

