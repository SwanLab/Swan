classdef LevelSetWithSphereInclusion < LevelSetCreator
    
     properties (Access = private)
        radius
        center        
    end
    
    methods (Access = public)
        
        function obj = LevelSetWithSphereInclusion(input)
            obj.compute(input);
        end
    end
    
    methods (Access = protected)
        
    function computeInitialLevelSet(obj)
        obj.computeRadius()
        obj.computeCircleCenter()
        obj.computeLevelSet()
        obj.computeDesignVariable()
    end
    
    end
    
    methods (Access = private)
        
        function computeRadius(obj)
            xLength = max(obj.nodeCoord.x) - min(obj.nodeCoord.x);
            yLength = max(obj.nodeCoord.y) - min(obj.nodeCoord.y);
            zLength = max(obj.nodeCoord.z) - min(obj.nodeCoord.z);
            maxInteriorRadius = min([xLength/2,yLength/2,zLength/2]);
            fracRadius = 1-1e-6;
            obj.radius = fracRadius*maxInteriorRadius;
        end
        
        function computeCircleCenter(obj)
            x = obj.nodeCoord.x;
            y = obj.nodeCoord.y;            
            z = obj.nodeCoord.z;
            
            obj.center(1) = 0.5*(max(x) + min(x));
            obj.center(2) = 0.5*(max(y) + min(y)); 
            obj.center(3) = 0.5*(max(z) + min(z));                         
        end
        
        function computeLevelSet(obj)
             x = obj.nodeCoord.x;
             y = obj.nodeCoord.y;
             z = obj.nodeCoord.z;             
             cx = obj.center(1);
             cy = obj.center(2);             
             cz = obj.center(3);                                       
             r = obj.radius;
             ls = ((x-cx)/r).^2 + ((y-cy)/r).^2 + ((z-cz)/r).^2 - 1;
             obj.levelSet = ls;
        end
        
        function computeDesignVariable(obj)
             switch obj.optimizerName
                case {'SLERP','HAMILTON-JACOBI'}
                    obj.x = obj.levelSet;
                otherwise
                    initial_holes = ceil(max(obj.levelSet,0))>0;
                    obj.x(initial_holes) = obj.hole_value;
            end
            
            
        end
        
    end
    
end


