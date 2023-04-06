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
            %R=0.25;
  
                x1  = obj.nodeCoord(:,1);
                x2  = obj.nodeCoord(:,2);
                mx1 = 1-x1.*x2;
                mx2 = 1-x1.*x2;
                xDef = obj.deformCoord(x1,x2);
               

                y1  = x1Def/obj.epsilon(1,1);
                y2  = x2Def/obj.epsilon(1,2);

                fpy1 = y1-floor(y1)-0.5;
                fpy2 = y2-floor(y2)-0.5;
                
                ls(:,1) = ((2*fpy1)./mx1).^q+((2*fpy2)./mx2).^q;

            obj.levelSet = 1-ls;
            
        end
        
        function init(obj,cParams)
            obj.epsilon = cParams.epsilon;
        end
        
    end
    
end
