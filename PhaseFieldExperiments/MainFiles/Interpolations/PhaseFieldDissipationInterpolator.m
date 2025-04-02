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

            s.operation = @(phi) abs(phi).^p;
            s.ndimf = 1;
            obj.fun = DomainFunction(s);

            s.operation = @(phi) p*(abs(phi)).^(p-1);
            s.ndimf = 1;
            obj.dfun = DomainFunction(s);

            s.operation = @(phi) p*(p-1)*(abs(phi+e)).^(p-2);
            s.ndimf = 1;
            obj.ddfun = DomainFunction(s);

        end

    end

end