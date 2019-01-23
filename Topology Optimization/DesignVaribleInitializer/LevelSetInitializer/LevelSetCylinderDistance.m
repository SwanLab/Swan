classdef LevelSetCylinderDistance < ...
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
            obj.computeDistanceToCylinder()
            obj.computeLevelSetValue()
        end
        
    end
    
    methods (Access = private)
                
        function computeDistanceToCylinder(obj)            
            r = obj.radius;
            d = zeros(obj.lsSize);
            for idim = 1:2
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

