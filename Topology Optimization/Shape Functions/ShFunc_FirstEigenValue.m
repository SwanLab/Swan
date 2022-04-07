classdef ShFunc_FirstEigenValue < ShapeFunctional
   
    properties (Access = private)

    end
    
    methods (Access = public)
        
        function obj = ShFunc_FirstEigenValue(cParams)
            obj.init(cParams)
        end
        
        function computeFunctionAndGradient(obj)
            obj.computeFunction();
            obj.computeGradient();
        end
        
    end

    methods (Access = protected)

        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
        end
        
    end
    
    methods (Access = public)

        function computeFunction(obj)
            gamma = obj.designVariable.getFirstEigenMode();
            f0val = -gamma;
            obj.value = f0val;
        end

        function computeGradient(obj)
            N = obj.designVariable.getNelem();
            df0dx = zeros(N+1,1);
            df0dx(N+1) = -1;
            obj.gradient = df0dx;   
        end

    end
    
end