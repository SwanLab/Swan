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

        function t = getTitlesToPlot(obj)
            t{1} = 'First Eigen Value';
        end            
        
        function v = getVariablesToPlot(obj)
            v{1} = obj.value;%*obj.value0;
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
            x = obj.designVariable.value;
            df0dx = zeros(size(x));
            df0dx(end) = -1;
            obj.gradient = df0dx;   
        end

    end
    
end