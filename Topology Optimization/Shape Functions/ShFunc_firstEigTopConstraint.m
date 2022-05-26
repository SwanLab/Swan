classdef ShFunc_firstEigTopConstraint < ShapeFunctional
    
    properties (Access = private)
        eigModesTopology
    end

    methods (Access = public)
        
        function obj = ShFunc_firstEigTopConstraint(cParams)
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
            v{1} = obj.value; % *obj.value0;
        end 
        
    end

    methods (Access = protected)

        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
            obj.eigModesTopology = cParams.eigModesTopology;
        end
    end
    
    methods (Access = public)

        function computeFunction(obj)
           obj.value = obj.eigModesTopology.provideFunction();
           %obj.normalizeFunction();
        end

        function computeGradient(obj)
            dfdx = obj.eigModesTopology.provideDerivative();
            obj.gradient = dfdx';
            %obj.normalizeGradient();
        end
    end
    
end