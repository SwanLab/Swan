classdef NumericalRectangleHomogenizer < NumericalRectangleTypeHomogenizer
    
    methods (Access = public)
        
        function obj = NumericalRectangleHomogenizer(fileName,print,m1,m2,iter)
            obj.compute(fileName,print,m1,m2,iter)            
        end
        
        
    end
 
    
    methods (Access = protected)
        
        function createLevelSet(obj,m1,m2)
            input.mesh = obj.microProblem.mesh;
            input.settings = obj.setting;
            input.epsilon = 0.1;
            input.m1 = m1;
            input.m2 = m2;
            input.coord = obj.microProblem.mesh;
            input.ndim = obj.microProblem.mesh.ndim;
            DesignVar = DesignVaribleInitializer_Rectangle(input); 
            obj.nodalLevelSet = DesignVar.x;
        end
        
    end
    
end

