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

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.val     = 1e-1;
            obj.oldVal  = 1e-1;
            obj.tolMin  = cParams.tolMin;
            obj.tolMax  = cParams.tolMax;
        end

        function compute(obj,indNorm)
            t       = indNorm^2;
            obj.val = max(obj.tolMin,max(t,obj.tolMin));
        end
        
        function updateOldTol(obj)
            obj.oldVal = obj.val;
        end

    end
    
end