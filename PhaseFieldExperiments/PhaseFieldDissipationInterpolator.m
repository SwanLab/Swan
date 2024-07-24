classdef PhaseFieldDissipationInterpolator < handle

    properties (Access = public)
        fun
        dfun
        ddfun
    end

    properties (Access = private)
        pExp
    end

    methods (Access = public)

        function obj = PhaseFieldDissipationInterpolator(cParams)
            obj.init(cParams)
            obj.computeDissipationFunctionAndDerivatives();
        end
        
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.pExp = cParams.pExp;
        end

        function computeDissipationFunctionAndDerivatives(obj)
            p = obj.pExp;
            e = 1e-12;

            obj.fun      = @(phi) abs(phi).^p;
            obj.dfun     = @(phi) p*(abs(phi)).^(p-1);
            obj.ddfun    = @(phi) p*(p-1)*(abs(phi+e)).^(p-2);
        end

    end

end