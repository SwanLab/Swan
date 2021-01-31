classdef ShFunc_Chomog_alphabeta < ShFunc_Chomog
    properties (Access = private)
        alpha
        beta
    end
    methods
        function obj=ShFunc_Chomog_alphabeta(cParams)
            obj@ShFunc_Chomog(cParams);
            obj.alpha = cParams.alpha/norm(cParams.alpha);
            obj.beta = cParams.beta/norm(cParams.beta);
        end
        function computeFunctionAndGradient(obj)
            obj.computePhysicalData();
            obj.computeFunctionValue();
            obj.normalizeFunction();            
            obj.computeGradient();
            obj.normalizeGradient();            
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