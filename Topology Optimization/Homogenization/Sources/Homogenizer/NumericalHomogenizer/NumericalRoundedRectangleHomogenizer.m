classdef NumericalRoundedRectangleHomogenizer < NumericalRectangleTypeHomogenizer
    
    
    methods (Access = public)
        
        function obj = NumericalRoundedRectangleHomogenizer(fileName,print,m1,m2,iter)
            obj.compute(fileName,print,m1,m2,iter)            
        end        
        
    end
 
    
    methods (Access = protected)
        
        function createLevelSet(obj,m1,m2)
            input.mesh = obj.microProblem.mesh;
            input.settings = obj.setting;
            input.epsilon = 0.1;
            DesignVar = DesignVaribleInitializerRoundedRectangle(input); 
            DesignVar.computeInitialLevelSet(m1,m2);
            obj.nodalLevelSet = DesignVar.x;
        end
        
    end
    
end

