classdef ShFunc_Chomog_alphabeta < ShFunc_Chomog
    
    properties (Access = private)
        alpha
        beta
    end
    
    methods (Access = public)
   
        function obj = ShFunc_Chomog_alphabeta(cParams)
            obj.initChomog(cParams);
            obj.alpha = cParams.alpha/norm(cParams.alpha);
            obj.beta  = cParams.beta/norm(cParams.beta);
        end
        
        function computeFunctionValue(obj)
            invCh   = inv(obj.Chomog);
            invChAB = obj.projectTensor(invCh,obj.alpha,obj.beta);
            obj.value = invChAB;
        end
        
        function computeGradientValue(obj)
            obj.computeChDerivative();
            dinvChAB = obj.computedChInv(obj.Chomog,obj.alpha,obj.beta);
            obj.gradient = dinvChAB;
        end
        
        function q = getQuad(obj)
            q = obj.physicalProblem.getQuadrature();
        end
        
    end
    

end