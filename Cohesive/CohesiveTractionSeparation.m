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
            t = obj.law.computeFunction(jump); %2x1 x(integration Point)
        end

        function d = computeDerivative(obj, jump)
            d = obj.law.computeDerivative(jump); %2x2 x(integration Point) matriu diagonal
        end

    end
    
    methods (Access = private)

    end
    
end