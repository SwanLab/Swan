classdef CohesiveCubicLaw < handle

    properties (Access = private)
        sigmaMax
        normalCharLength
        tangencialCharLength
    end
    
    methods (Access = public)
        
        function obj = CubicCohesiveLaw(cParams)
            obj.init(cParams)            
        end

        function t = computeFunction(obj,disp)
            % norm(disp)
            t = 6.75*obj.sigmaMax*disp*(1-2*disp+disp^2); % no se si aixo esta be
        end

        function d = computeDerivative(obj,disp)
            d = 6.75 * obj.sigmaMax * (1 - 4*disp + 3*disp^2);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.sigmaMax = cParams.sigmaMax;
            obj.normalCharLength = cParams.normalCharLength;
            obj.tangencialCharLength = cParams.tangencialCharLength;
        end
        
    end
    
end