classdef ShFunc_FirstEigenValue < ShapeFunctional
   
    properties (Access = private)
        nElem
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
            obj.nElem = cParams.nElem;
            obj.designVariable = cParams.designVariable;
        end
        
    end
    
    methods (Access = public)

        function computeFunction(obj)
            N = obj.nElem;
            x = obj.designVariable.value;
            f0val = -x(N+1);
            obj.value = f0val;
        end

        function computeGradient(obj)
            N = obj.nElem;
            df0dx = zeros(N+1,1);
            df0dx(N+1) = -1;
            obj.gradient = df0dx;   
        end

    end
    
end