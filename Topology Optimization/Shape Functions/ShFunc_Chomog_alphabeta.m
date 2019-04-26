classdef ShFunc_Chomog_alphabeta < ShFunc_Chomog
    properties (Access = private)
        alpha
        beta
    end
    methods
        function obj=ShFunc_Chomog_alphabeta(settings)
            obj@ShFunc_Chomog(settings);
            obj.alpha=settings.alpha/norm(settings.alpha);
            obj.beta=settings.beta/norm(settings.beta);
        end
        function computeCostAndGradient(obj)
            obj.computePhysicalData();
            obj.computeFunctionValue();
            obj.computeGradient();
            obj.normalizeFunctionAndGradient();
        end
    end
    
    methods (Access = private)
        
        function computeGradient(obj)
            obj.compute_Chomog_Derivatives();
            inv_matCh = inv(obj.Chomog);
            gradient = obj.derivative_projection_Chomog(inv_matCh,obj.alpha,obj.beta);            
            mass     = obj.Msmooth;
            gradient = obj.filter.getP1fromP0(gradient');
            gradient = mass*gradient;
            obj.gradient = gradient;
        end
        
        function computeFunctionValue(obj)
            inv_matCh = inv(obj.Chomog);
            c = obj.projection_Chomog(inv_matCh,obj.alpha,obj.beta);
            obj.value = c;
        end
        
    end
end