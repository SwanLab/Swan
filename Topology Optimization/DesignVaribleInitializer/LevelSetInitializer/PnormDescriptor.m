classdef PnormDescriptor < NormDescriptor
    
    properties (Access = protected)
        pnorm
    end
    
    methods (Access = protected)
        
        function d = computeDistance(obj)
            x = obj.pos;
            p = obj.pnorm;
            xnorm = zeros(size(x,1),1);
            for idim = 1:size(x,2)
               xnorm = xnorm + x(:,idim).^p; 
            end
            d = xnorm.^(1/p) + 1e-14 - 1;
            obj.dist = d;
        end
        
    end
    
end

