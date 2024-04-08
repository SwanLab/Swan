classdef ConjugateGradientToleranceCalculator < handle

    properties (Access = public)
        val
        oldVal
    end

    properties (Access = private)
        tolMax,tolMin
    end

    methods (Access = public)

        function obj = ConjugateGradientToleranceCalculator(cParams)
            obj.init(cParams);
        end

        function compute(obj,indNorm)
            t          = indNorm;
            t          = 0.5*(t + obj.oldVal); % Relaxing the change of TOL
            obj.val    = max(obj.tolMin,min(t,obj.tolMax));
            obj.oldVal = obj.val;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.val     = 1e-3;
            obj.oldVal  = 1e-3;
            obj.tolMin  = cParams.tolMin;
            obj.tolMax  = cParams.tolMax;
        end
        
        function updateOldTol(obj)
            obj.oldVal = obj.val;
        end

    end

end