classdef PnormDescriptor < NormDescriptor
    
    properties (Access = protected)
        pnorm
    end
    
    properties (Access = protected, Abstract)
       pos 
    end
    
    methods (Access = protected)
        
        function computeDistance(obj)
            x = obj.pos;
            p = obj.pnorm;
            xnorm = zeros(size(x,1),1);
            for idim = 1:size(x,2)
               xnorm = xnorm + abs(x(:,idim)).^p; 
            end
            d = xnorm.^(1./p);
            obj.dist = d;
        end
        
    end
    
end

