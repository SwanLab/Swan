classdef PhaseFieldIsotropicHomogenizationDegradation < handle

    properties (Access = private)
        shearFun
        bulkFun
    end

    properties (Access = private)
        fileName
        E
        mesh
    end

    methods (Access = public)
        function obj = PhaseFieldIsotropicHomogenizationDegradation(cParams)
            obj.init(cParams);
            obj.computeShearAndBulkFromData();
        end

        function [mu,kappa] = computeConstitutiveTensor(obj,phi)
            mu    = obj.interpolate(obj.shearFun,phi);
            kappa = obj.interpolate(obj.bulkFun,phi);
        end

        function [dmu,dkappa] = computeConstitutiveTensorDerivative(obj,phi)
            dmu    = obj.derive(obj.shearFun,phi);
            dkappa = obj.derive(obj.bulkFun,phi);
        end

        function [ddmu,ddkappa] = computeConstitutiveTensorSecondDerivative(obj,phi)
            ddmu    = obj.derive2(obj.shearFun,phi);
            ddkappa = obj.derive2(obj.bulkFun,phi);
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fileName = cParams.params.fileName;
            obj.E        = cParams.young.constant;
            obj.mesh     = cParams.young.mesh;
        end

        function computeShearAndBulkFromData(obj)
            [data] = load(obj.fileName);
            C11 = data.degradation.fun{1,1,1,1};
            C33 = data.degradation.fun{1,2,1,2};
            obj.shearFun = @(phi) obj.E.*(C11(phi)-C33(phi));
            obj.bulkFun  = @(phi) obj.E.*C33(phi);
        end

        function f = interpolate(obj,fun,phi)
            f = obj.createDomainFunction(fun,phi);
        end

        function df = derive(obj,fun,phi)
            syms x
            dfun = matlabFunction(diff(fun(x)));
            df = obj.createDomainFunction(dfun,phi);
        end

        function d2f = derive2(obj,fun,phi)
            syms x
            d2fun = matlabFunction(diff(diff(fun(x))));
            d2f = obj.createDomainFunction(d2fun,phi);
        end

        function f = createDomainFunction(obj,fun,phi)
            s.operation = @(xV) fun(phi.evaluate(xV));
            s.ndimf = 1;
            s.mesh = obj.mesh;
            f = DomainFunction(s);
        end

    end

end