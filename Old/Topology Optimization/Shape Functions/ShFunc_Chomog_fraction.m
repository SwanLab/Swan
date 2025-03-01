classdef ShFunc_Chomog_fraction < ShFunc_Chomog
   
    properties (Access = private)
        invChAA
        invChBB
        invChAB
        invChBA
    end
    
    properties (Access = private)
        alpha
        beta
    end
    
    methods (Access = public)
        
        function obj=ShFunc_Chomog_fraction(cParams)
            obj.initChomog(cParams);
            obj.alpha = cParams.alpha/norm(cParams.alpha);
            obj.beta  = cParams.beta/norm(cParams.beta);
        end
        
        function computeFunctionValue(obj) 
            obj.computeInvChProyections();
            obj.value = obj.invChAB/obj.invChAA + obj.invChBA/obj.invChBB;
        end
        
        function computeGradientValue(obj)
            obj.computeChDerivative();
            a = obj.alpha;
            b = obj.beta;
            beta1 = obj.invChAA*b - obj.invChAB*a;
            beta2 = obj.invChBB*a - obj.invChBA*b;
            g1    = obj.computedChInv(obj.Chomog,a,beta1);
            g2    = obj.computedChInv(obj.Chomog,b,beta2);
            grad = g1/(obj.invChAA)^2 + g2/(obj.invChBB)^2;
            obj.gradient = grad;
        end
        
    end
    
    methods (Access = private)
        
        function computeInvChProyections(obj)
            invCh = inv(obj.Chomog);
            a = obj.alpha;
            b = obj.beta;
            obj.invChAB  = obj.projectTensor(invCh,a,b);
            obj.invChBA  = obj.projectTensor(invCh,b,a);
            obj.invChAA  = obj.projectTensor(invCh,a,a);
            obj.invChBB  = obj.projectTensor(invCh,b,b);
        end
        
    end

end