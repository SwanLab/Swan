classdef CohesiveCubicLaw < handle

    properties (Access = private)
        sigmaMax
        normalCharLength
        tangencialCharLength
    end
    
    methods (Access = public)
        
        function obj = CohesiveCubicLaw(cParams)
            obj.init(cParams)            
        end

        function t = computeFunction(obj,jump)

            t = 6.75.*obj.sigmaMax.*jump.*(1-2.*jump+jump.^2);
            
            
            % t es una domain 2x1
        end

        function d = computeDerivative(obj,jump)


            d = 6.75 * obj.sigmaMax * (1 - 4*jump + 3*jump^2);


        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.sigmaMax                = cParams.sigmaMax;
            obj.normalCharLength        = cParams.normalCharLength;
            obj.tangencialCharLength    = cParams.tangencialCharLength;
        end

        function eff = computeEffectiveJumps(obj,jump)
            eff = norm( (1/obj.tangencialCharLength;0) * jump)


        end
        
    end
    
end