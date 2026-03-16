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

        function t = evaluate(obj,xV,jump)

            t = 6.75.*obj.sigmaMax.*jump(xV).*(1-2.*jump(xV)+jump(xV).^2);

        end

        function d = derivative(obj,jump)
            d = 6.75 * obj.sigmaMax * (1 - 4*jump + 3*jump^2);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.sigmaMax                = cParams.sigmaMax;
            obj.normalCharLength        = cParams.normalCharLength;
            obj.tangencialCharLength    = cParams.tangencialCharLength;
        end
        
    end
    
end