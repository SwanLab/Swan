classdef LevelSetYounes < LevelSetCreator

    properties (Access = private)
        epsilon
    end
    
    methods (Access = public)
        
        function obj = LevelSetYounes(cParams)
            obj.init(cParams);
            obj.compute(cParams);
        end

        function mF = getFineMesh(obj)
            fMesh = obj.remesher.fineMesh;
            mF = fMesh.createDiscontinuousMesh();
        end
        
    end
    
    methods (Access = protected)
        
        function computeLevelSet(obj)
            obj.computeLevelSetValue();
        end
        
    end
    
    methods (Access = private)
        
        function computeLevelSetValue(obj)
            q=2;
            for i=1:size(obj.nodeCoord,1)
                x1  = obj.nodeCoord(i,1);
                x2  = obj.nodeCoord(i,2);
                y1  = x1/obj.epsilon(1,1);
                y2  = x2/obj.epsilon(1,2);
                fpy1 = y1-floor(y1)-0.5;
                fpy2 = y2-floor(y2)-0.5;
                mx1 = 1-x1*x2;
                mx2 = 1-x1*x2;

                ls(i,1) = ((2*fpy1)/mx1)^q+((2*fpy2)/mx2)^q-1;
            end
            obj.levelSet = ls;
            
        end
        
        function init(obj,cParams)
            obj.epsilon = cParams.epsilon;
        end
        
    end
    
end
