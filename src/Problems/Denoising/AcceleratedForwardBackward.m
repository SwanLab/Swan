classdef AcceleratedForwardBackward < SplittingAlgorithm

    properties (Access = private)
        gradientMethod
        proximal
        momentum
    end
    
    methods (Access = public)
        
        function obj = AcceleratedForwardBackward(cParams)
            obj.proximal = cParams.proximal;
            obj.createGradientMethods(cParams);
            obj.createMomentum(cParams);
        end
        
        function update(obj)
            obj.momentum.apply();
            obj.gradientMethod.compute();
            obj.proximal.solve();
        end
        
    end
    
    methods (Access = private)
        
        function createGradientMethods(obj,cParams)
            s = cParams.gradientMethodParams;
            obj.gradientMethod = GradientMethod(s);
        end
        
        function createMomentum(obj,cParams)
            s = cParams.momentumParams;
            obj.momentum = Momentum(s);
        end
       
    end
    
end