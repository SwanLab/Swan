classdef Sh_doubleEig < ShapeFunctional
    
    properties (Access = private)
        eigNum
        eigModes
        name
    end

    methods (Access = public)
        
        function obj = Sh_doubleEig(cParams)
            obj.init(cParams)
        end

        function computeFunctionAndGradient(obj)
            obj.computeFunction();
            obj.computeGradient();
        end

        function t = getTitlesToPlot(obj)
            t{1} = obj.name;
        end  

        function v = getVariablesToPlot(obj)
            v{1} = obj.value; % *obj.value0;
        end 
        
    end

    methods (Access = protected)

        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
            obj.eigModes = cParams.eigModes;
            obj.eigNum = cParams.eigNum;
            obj.name = cParams.name;
        end
    end
    
    methods (Access = public)

        function computeFunction(obj)
            eigN = obj.eigNum;
            lambda = obj.eigModes.provideEigenValue();
            gamma = obj.designVariable.getFirstEigenMode();            
            obj.value = gamma-lambda(eigN);
        end

        function computeGradient(obj)
            eigN = obj.eigNum;
            dfdx = obj.eigModes.provideDerivative(eigN);
            obj.gradient = dfdx';
        end
    end
    
end