classdef ShFunc_Perimeter < Shape_Functional
    
    properties (Access = protected)
        epsilon
        designVariable
        regularizedDensity
        regularizedDensityProjection
    end
    
    methods
        function obj = ShFunc_Perimeter(settings)
            settings.filterParams.filterType = 'PDE';
            obj.init(settings);
        end
        
        function computeCostAndGradient(obj,designVariable)
            obj.updateProtectedVariables(designVariable)
            obj.computeRegularizedDensity()
            obj.computeRegularizedDensityProjection()
            obj.computePerimeterValue()
            obj.computePerimeterGradient()
        end
                
    end
    
    methods (Access = protected)
        
        function updateProtectedVariables(obj,designVariable)
            obj.updateDesignVariable(designVariable)
            obj.updateEpsilonValue()
            obj.updateEpsilonInFilter()
        end
        
        function updateDesignVariable(obj,designVariable)
            obj.designVariable = designVariable;
        end
       
        function updateEpsilonValue(obj)
            obj.epsilon=obj.target_parameters.epsilon_perimeter;
        end
        
        function updateEpsilonInFilter(obj)
            obj.filter.updateEpsilon(obj.epsilon);
        end
        
        function computeRegularizedDensity(obj)
            obj.regularizedDensity = obj.filter.getP1fromP1(obj.designVariable);
        end
        
        function computeRegularizedDensityProjection(obj)
            obj.regularizedDensityProjection = obj.filter.integrate_L2_function_with_shape_function(obj.designVariable);
        end
        
        function computePerimeterValue(obj)
            obj.value = 0.5/obj.epsilon*((1 - obj.regularizedDensity)'*obj.regularizedDensityProjection);
        end
        
        function computePerimeterGradient(obj)
            obj.computeContinousGradient();
            obj.computeDiscreteGradient();
        end
        
        function computeContinousGradient(obj)
            obj.gradient = 0.5/obj.epsilon*(1 - 2*obj.regularizedDensity);
        end
        
        function computeDiscreteGradient(obj)
            obj.gradient = obj.Msmooth*obj.gradient;
        end
        
    end
end