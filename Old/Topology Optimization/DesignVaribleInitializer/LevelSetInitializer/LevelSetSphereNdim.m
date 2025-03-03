classdef LevelSetSphereNdim < ...
        LevelSetCreator & ...
        LevelSetCenterDescriptor & ...
        LevelSetRadiusDescriptor
    
    properties (Access = protected)  
        fracRadius
        dist
    end
    
    methods (Access = protected)
        
        function computeLevelSet(obj)
            obj.computeRadius()
            obj.computeCenter()
            obj.computeDistanceToSphere()
            obj.computeLevelSetValue()
        end
        
    end
    
    methods (Access = private)
                
        function computeDistanceToSphere(obj)            
            r = obj.radius;
            d = zeros(obj.lsSize);
            for idim = 1:obj.ndim
                pos = obj.nodeCoord(:,idim);
                pos0 = obj.center(idim);
                d = d + ((pos-pos0)/r).^2;
            end
            obj.dist = d;
        end
        
    end

    methods (Access = protected)
       computeLevelSetValue(obj) 
    end
    
end

