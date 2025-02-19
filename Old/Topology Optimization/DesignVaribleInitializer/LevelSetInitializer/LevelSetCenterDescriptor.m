classdef LevelSetCenterDescriptor < handle
    
    properties (Access = protected)
        center
    end
    
    properties (Access = protected, Abstract)
        nodeCoord
        ndim
    end
    
    methods (Access = protected)
        
        function computeCenter(obj)
            for idim = 1:obj.ndim
                pos = obj.nodeCoord(:,idim);
                obj.center(idim) = 0.5*(max(pos) + min(pos));
            end
        end
        
    end
    
end

