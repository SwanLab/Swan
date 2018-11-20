classdef LevelSetSphereNdim < LevelSetCreator
    
     properties (Access = protected)
        radius
        center  
        fracRadius
    end
    
    methods (Access = protected)
        
    function computeInitialLevelSet(obj)
        obj.computeRadius()
        obj.computeCenter()
        obj.computeLevelSet()
        obj.computeDesignVariable()
    end
    
    end
    
    methods (Access = private)
        
        function computeRadius(obj)
            lengthDim = zeros(obj.ndim,1);
            for idim = 1:obj.ndim
                pos = obj.nodeCoord(:,idim);
                lengthDim(idim) = 0.5*(max(pos) - min(pos));
            end
            maxInteriorRadius = min(lengthDim);            
            obj.radius = obj.fracRadius*maxInteriorRadius;
        end
        
        function computeCenter(obj)
            for idim = 1:obj.ndim
                pos = obj.nodeCoord(:,idim);
                obj.center(idim) = 0.5*(max(pos) + min(pos));
            end
        end
        
        function computeLevelSet(obj)
            r = obj.radius;
            ls = zeros(obj.lsSize);
            for idim = 1:obj.ndim
                pos = obj.nodeCoord(:,idim);
                pos0 = obj.center(idim);
                ls = ls + ((pos-pos0)/r).^2;
            end
            ls = ls - 1;
            obj.levelSet = ls;
        end
        
    end
    
    methods (Abstract, Access = protected)
        computeDesignVariable(obj)
    end
    
end

