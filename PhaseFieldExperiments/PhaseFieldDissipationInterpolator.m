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
            
            obj.fun      = @(phi) phi.^p;
            obj.dfun     = @(phi) p*(phi).^(p-1);
            if p == 1
                obj.ddfun  = @(phi) 0*(phi);
            else
                obj.ddfun  = @(phi) p*(p-1)*(phi).^(p-2);
            end
        end

    end

end