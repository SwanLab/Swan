classdef MaterialPhaseField < IsotropicElasticMaterial
    
    properties (Access = public)
        ndimf
    end

    properties (Access = private)
        fun
        phi
        tensorType
    end

    properties (Access = private)
        mesh
        matInterpolation
        Gc
    end

    methods (Access = public)

        function obj = setMaterial(obj,phi,interpType,tensorType)
            obj.phi = phi;
            obj.tensorType = tensorType;
            switch interpType
                case 'Isotropic'
                    obj.fun = @(phi) 1;
                case 'Interpolated'
                    obj.fun = obj.matInterpolation.fun;
                case 'Jacobian'
                    obj.fun = obj.matInterpolation.dfun;
                case 'Hessian'
                    obj.fun = obj.matInterpolation.ddfun;
            end
        end

        function obj = MaterialPhaseField(cParams)
            obj.init(cParams)
            obj.initPhaseField(cParams)
        end

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

        function kFun = getBulkFun(obj,phi,interpType)
            obj.setMaterial(phi,interpType,'');
            g = obj.computeDegradationFun();

            E = obj.young.fValues;
            nu = obj.poisson.fValues;
            kV = obj.computeKappaFromYoungAndPoisson(E,nu,obj.ndim);

            s.mesh = obj.mesh;
            s.order = obj.young.order;
            s.fValues = kV;
            k = LagrangianFunction(s);

            kFun = times(g,k);
        end

        function muFun = getShearFun(obj,phi,interpType)
            obj.setMaterial(phi,interpType,'');
            g = obj.computeDegradationFun();

            E = obj.young.fValues;
            nu = obj.poisson.fValues;
            muV = obj.computeMuFromYoungAndPoisson(E,nu);
            
            s.mesh = obj.mesh;
            s.order = 'P1';
            s.fValues = muV;
            mu = LagrangianFunction(s);

            muFun = g.*mu;
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

        function [mu,l] = computeMatTensorParams(obj,xV)
            g = obj.computeDegradationFun();
            [mu, k] = obj.computeMuAndKappa(g,xV);
            l = obj.computeLambdaFromShearAndBulk(mu,k,obj.ndim);
        end

        function g = computeDegradationFun(obj)
            s.operation = @(xV) obj.fun(obj.phi.evaluate(xV));
            s.ndimf = obj.phi.ndimf;
            g = DomainFunction(s);
        end

        function [mu, k] = computeMuAndKappa(obj,g,xV)
            [muV,kV] = obj.computeShearAndBulk(xV);
            if obj.tensorType == "Deviatoric"
                kV = 0.*kV;
            elseif obj.tensorType == "Volumetric"
                muV = 0.*muV; 
            end

            mu = g.evaluate(xV).*muV;
            k = g.evaluate(xV).*kV;
        end
    end

end