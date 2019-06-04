classdef ShFunc_Perimeter < ShapeFunctional
    
    properties (Access = protected)
        epsilon
        regularizedDensity
        regularizedDensityProjection
    end
    
    methods
        function obj = ShFunc_Perimeter(cParams)
            cParams.filterParams.quadratureOrder = 'LINEAR';            
            cParams.filterParams.filterType = 'PDE';
            obj.init(cParams);
        end
        
        function computeCostAndGradient(obj)
            obj.updateProtectedVariables()
            obj.computeRegularizedDensity()
            obj.computeRegularizedDensityProjection()
            obj.computePerimeterValue()
            obj.computePerimeterGradient()
        end
                
    end
    
    methods (Access = protected)
        
        function updateProtectedVariables(obj)
            obj.updateEpsilonValue()
            obj.updateEpsilonInFilter()
        end        
       
        function updateEpsilonValue(obj)
            obj.epsilon = obj.target_parameters.epsilon_perimeter;
        end
        
        function updateEpsilonInFilter(obj)
            obj.filter.updateEpsilon(obj.epsilon);
        end
        
        function computeRegularizedDensity(obj)
            obj.regularizedDensity = obj.filter.getP1fromP1(obj.designVariable.value);
        end
        
        function computeRegularizedDensityProjection(obj)
            obj.regularizedDensityProjection = obj.filter.integrate_L2_function_with_shape_function(obj.designVariable.value);
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