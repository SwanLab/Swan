classdef Sh_firstEigTop < ShapeFunctional
    
    methods (Access = public)
        
        function obj = Sh_firstEigTop(cParams)
            obj.init(cParams)
        end
        
        function computeFunctionAndGradient(obj)
            obj.computeFunctionValue();
            obj.computeGradientValue();
        end

        function t = getTitlesToPlot(obj)
            t{1} = 'First Eigen Value';
        end            
        
        function v = getVariablesToPlot(obj)
            v{1} = obj.value*obj.value0;
        end
        
    end

    methods (Access = protected)
        
        function computeFunctionValue(obj) % to discuss the formulation
            obj.value
            
        end
        
        function computeGradientValue(obj) % to discuss the formulation
            obj.gradient
        end
        
    end
    
end