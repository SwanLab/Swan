classdef CohesiveTractionSeparation < handle
    
    properties (Access = public)
        lawType
        law
    end
    
    methods (Access = public)
        
        function obj = CohesiveTractionSeparation(cParams)
            obj.law = CohesiveLawFactory.create(cParams); 
        end

        function t = evaluate(obj, xV, jump)
            t = obj.law.evaluate(xV,jump); %2x1 x(integration Point)
        end

        function d = derivative(obj,xV, jump)
            d = obj.law.derivative(xV,jump); %2x2 x(integration Point) matriu diagonal
        end

        
    end
    
    methods (Access = private)

    end
    
end