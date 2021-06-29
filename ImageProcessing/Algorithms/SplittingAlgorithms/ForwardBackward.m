classdef ForwardBackward < SplittingAlgorithm
    
    properties (Access = private)
        gradientMethod
        proximal
    end
    
    methods (Access = public)
        
        function obj = ForwardBackward(cParams)
            obj.proximal = cParams.proximal;
            obj.createGradientMethods(cParams);
        end
        
        function update(obj)
            obj.gradientMethod.compute();
            obj.proximal.solve();
        end
        
    end
    
    methods (Access = private)
        
        function createGradientMethods(obj,cParams)
            s = cParams.gradientMethodParams;
            obj.gradientMethod = GradientMethod(s);
        end
       
    end
    
    
end