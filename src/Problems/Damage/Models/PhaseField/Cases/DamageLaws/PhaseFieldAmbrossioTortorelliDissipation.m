classdef PhaseFieldAmbrossioTortorelliDissipation < handle

    properties (Access = public)
        interpolation
        constant
    end

    properties (Access = private)
        mesh
        pExp
    end

    methods (Access = public)

        function obj = PhaseFieldAmbrossioTortorelliDissipation(cParams)
            obj.init(cParams)
            obj.computeDissipationFunctionAndDerivatives();
            obj.computeDissipationConstant(cParams);
        end
        
    end

    methods (Access = private)

        function init(obj,cParams)
            if ~ismember(cParams.pExp, [1 2])
                error('Exponent must be 1 or 2');
            else
                obj.pExp = cParams.pExp;
                obj.mesh = cParams.mesh;
            end
        end

        function computeDissipationConstant(obj,cParams)
            if isfield(cParams,'constant')
                obj.constant = cParams.constant;
            else
                obj.constant = 4*(2/(obj.pExp+2));
            end
        end

        function computeDissipationFunctionAndDerivatives(obj)
            p = obj.pExp;
            e = 1e-12;

            s.operation = @(phi) abs(phi+e).^p;
            s.ndimf = 1;
            s.mesh = obj.mesh;
            obj.interpolation.fun = DomainFunction(s);

            s.operation = @(phi) p*(abs(phi+e)).^(p-1);
            s.ndimf = 1;
            s.mesh = obj.mesh;
            obj.interpolation.dfun = DomainFunction(s);

            s.operation = @(phi) p*(p-1)*(abs(phi+e)).^(p-2);
            s.ndimf = 1;
            s.mesh = obj.mesh;
            obj.interpolation.ddfun = DomainFunction(s);
        end

    end

end