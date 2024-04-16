classdef MaterialPhaseField < IsotropicElasticMaterial
    
    properties (Access = public)
        ndimf
    end

    properties (Access = private)
        fun
        phi
        u
        tensorType
    end

    properties (Access = private)
        mesh
        matInterpolation
        Gc
    end

    methods (Access = public)

        function obj = MaterialPhaseField(cParams)
            obj.init(cParams)
            obj.initPhaseField(cParams)
        end

        function C = obtainTensor(obj)
            obj.fun = obj.matInterpolation.fun;
            s.operation = @(xV) obj.evaluate(xV);
            s.ndimf = 9;
            C{1} = DomainFunction(s);
        end

        function dC = obtainTensorDerivative(obj)
            obj.fun = obj.matInterpolation.dfun;
            s.operation = @(xV) obj.evaluate(xV);
            s.ndimf = 9;
            dC{1} =  DomainFunction(s);

        end

        function ddC = obtainTensorSecondDerivative(obj)
            obj.fun = obj.matInterpolation.ddfun;
            s.operation = @(xV) obj.evaluate(xV);
            s.ndimf = 9;
            ddC{1} =  DomainFunction(s);
        end

        function kFun = getBulkFun(obj,u,phi,interpType)
            obj.setMaterial(u,phi,interpType,'');
            [g,~] = obj.computeDegradationFun();

            E = obj.young.fValues;
            nu = obj.poisson.fValues;
            kV = obj.computeKappaFromYoungAndPoisson(E,nu,obj.ndim);

            s.mesh = obj.mesh;
            s.order = obj.young.order;
            s.fValues = kV;
            k = LagrangianFunction(s);

            kFun = times(g,k);
        end

        function muFun = getShearFun(obj,u,phi,interpType)
            obj.setMaterial(u,phi,interpType,'');
            [~,g0] = obj.computeDegradationFun();

            E = obj.young.fValues;
            nu = obj.poisson.fValues;
            muV = obj.computeMuFromYoungAndPoisson(E,nu);
            
            s.mesh = obj.mesh;
            s.order = 'P1';
            s.fValues = muV;
            mu = LagrangianFunction(s);

            muFun = g0.*mu;
        end

        function obj = setDesignVariable(obj,u,phi,tensorType)
            obj.u = u;
            obj.phi = phi;
            obj.tensorType = tensorType;
        end

    end

    methods (Access = private)

        function initPhaseField(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.matInterpolation = cParams.materialInterpolation;
            obj.Gc = cParams.Gc;
            obj.ndimf = 9;
        end

    end

    methods (Access = private)

        function C = evaluate(obj,xV)
            [mu,l] = obj.computeMatTensorParams(xV);
            nPoints = size(mu,3);
            nElem = size(mu,4);
            nStre = 3;
            C = zeros(nStre,nStre,nPoints,nElem);
            C(1,1,:,:)= 2*mu+l;
            C(1,2,:,:)= l;
            C(2,1,:,:)= l;
            C(2,2,:,:)= 2*mu+l;
            C(3,3,:,:)= mu;
        end

        function [mu,l] = computeMatTensorParams(obj,xV)
            [g, g0] = obj.computeDegradationFun();
            [mu, k] = obj.computeMuAndKappa(g,g0,xV);
            l = obj.computeLambdaFromShearAndBulk(mu,k,obj.ndim);
        end

        function [g, g0] = computeDegradationFun(obj)
            s.operation = @(xV) obj.fun(obj.phi.evaluate(xV));
            s.ndimf = obj.phi.ndimf;
            g0 = DomainFunction(s);

            trcSign = Heaviside(trace(SymGrad(obj.u)));
            g = g0.*trcSign + (1-trcSign);
        end

        function [mu, k] = computeMuAndKappa(obj,g,g0,xV)
            [muV,kV] = obj.computeShearAndBulk(xV);
            if obj.tensorType == "Deviatoric"
                kV = 0.*kV;
            elseif obj.tensorType == "Volumetric"
                muV = 0.*muV; 
            end
            mu = g0.evaluate(xV).*muV;
            k  = g.evaluate(xV).*kV;
        end
    end

end