classdef NumericalRectangleHomogenizer < NumericalRectangleTypeHomogenizer
    
    methods (Access = public)
        
        function obj = NumericalRectangleHomogenizer(fileName,print,m1,m2,iter)
            obj.compute(fileName,print,m1,m2,iter)            
        end
        
        
    end
 
    
    methods (Access = protected)
        
        function createLevelSet(obj,m1,m2)
            mesh = obj.microProblem.mesh;
            DesignVar = DesignVaribleInitializer_Rectangle(obj.setting,mesh,0.1); 
            DesignVar.compute_initial_x(m1,m2);
            obj.nodalLevelSet = DesignVar.x;
        end
        
    end
    
end

