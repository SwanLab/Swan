classdef LevelSetRadiusDescriptor < handle
    
    properties (Access = protected)
        radius
    end
    
    properties (Access = protected, Abstract)
        nodeCoord
        ndim
        fracRadius
    end
    
    methods (Access = protected)
        
        function computeRadius(obj)
            lengthDim = zeros(obj.ndim,1);
            for idim = 1:obj.ndim
                pos = obj.nodeCoord(:,idim);
                lengthDim(idim) = 0.5*(max(pos) - min(pos));
            end
            maxInteriorRadius = min(lengthDim);
            obj.radius = obj.fracRadius*maxInteriorRadius;
        end
        
    end
    
end

