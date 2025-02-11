classdef Sh_doubleEig < ShapeFunctional
    
    properties (Access = private)
        eigNum
        eigModes
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
            t{1} = 'Double Eigen Value';
        end  

        function v = getVariablesToPlot(obj)
            v{1} = obj.value; %*obj.value0;
        end 
        
    end

    methods (Access = protected)

        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
            obj.eigModes = cParams.eigModes;
            obj.eigNum = cParams.eigNum;
        end
    end
    
    methods (Access = public)

        function computeFunction(obj)
           eigN = obj.eigNum;
           obj.value = obj.eigModes.provideFunction(eigN);
           %obj.normalizeFunction();
        end

        function computeGradient(obj)
            eigN = obj.eigNum;
            dfdx = obj.eigModes.provideDerivative(eigN);
            obj.gradient = dfdx';
            %obj.normalizeGradient();
        end
    end
    
end