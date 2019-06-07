classdef ShFunc_Chomog_fraction < ShFunc_Chomog
   
    properties (Access = private)
        alpha
        beta
        
        Saa
        Sbb
        Sab
        Sba
        S
    end
    
    methods (Access = public)
        
        function obj=ShFunc_Chomog_fraction(settings)
            obj@ShFunc_Chomog(settings);
            obj.alpha=settings.alpha/norm(settings.alpha);
            obj.beta=settings.beta/norm(settings.beta);
        end
        
        function computeCostAndGradient(obj)
            obj.computePhysicalData();
            obj.computeInvChProyections();
            obj.computeFunctionValue();
            obj.computeGradient();
            obj.normalizeFunctionAndGradient();
        end
    end
    
    methods (Access = private)
                
        function computeGradient(obj)
            obj.compute_Chomog_Derivatives(); 
            a = obj.alpha;
            b = obj.beta;
            beta1 = obj.Saa*b - obj.Sab*a;
            beta2 = obj.Sbb*a - obj.Sba*b;
            g1 = obj.derivative_projection_Chomog(obj.S,a,beta1);
            g2 = obj.derivative_projection_Chomog(obj.S,b,beta2);
            gradient = g1/(obj.Saa)^2 + g2/(obj.Sbb)^2;            
            mass = obj.Msmooth;
            gradient=obj.filter.getP1fromP0(gradient(:));
            gradient = mass*gradient;
            obj.gradient = gradient;
        end
        
        function computeFunctionValue(obj)            
            obj.value = obj.Sab/obj.Saa + obj.Sba/obj.Sbb;
        end   
        
        function computeInvChProyections(obj)
            obj.S = inv(obj.Chomog);
            a = obj.alpha;
            b = obj.beta;
            obj.Sab  = obj.projection_Chomog(obj.S,a,b);
            obj.Sba  = obj.projection_Chomog(obj.S,b,a);
            obj.Saa  = obj.projection_Chomog(obj.S,a,a);
            obj.Sbb  = obj.projection_Chomog(obj.S,b,b);
        end
        
    end

end