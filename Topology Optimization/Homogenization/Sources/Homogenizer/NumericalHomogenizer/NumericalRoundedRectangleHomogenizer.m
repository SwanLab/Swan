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
            input.m1 = m1;
            input.m2 = m2;
            input.coord = obj.microProblem.mesh.coord;
            input.ndim  = obj.microProblem.mesh.ndim;
            designVar = LevelSetRectangleInclusion(input); 
            obj.nodalLevelSet = designVar.x;
        end
        
    end
    
end

