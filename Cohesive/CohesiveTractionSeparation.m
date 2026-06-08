classdef CohesiveTractionSeparation < handle
    
    properties (Access = public)
        lawType
        law
    end
    
    methods (Access = public)
        
        function obj = CohesiveTractionSeparation(cParams)
            obj.law = CohesiveLawFactory.create(cParams); 
        end

        function t = computeFunction(obj, jump)
            t = obj.law.computeFunction(jump);
        end

        function d = computeDerivativeTangent(obj, jump)
            d = obj.law.computeDerivativeTangent(jump); 
        end

        function d = computeDerivativeSecant(obj, jump)
            d = obj.law.computeDerivativeSecant(jump);
        end

        function dmgValues = getDamageValues(obj,jump)
            dmgValues = obj.law.computeDamage(jump);
        end

        function updateDamageOld(obj,jump)
            obj.law.updateDamageOld(jump);
        end

    end
    
    methods (Access = private)

    end
    
end