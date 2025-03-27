classdef MaterialPhaseFieldAnalyticalSplit < Material

    properties (Access = private)
        materialInterpolator
        mesh
    end

    methods (Access = public)

        function obj = MaterialPhaseFieldAnalyticalSplit(cParams)
            obj.init(cParams)
        end


        function V = obtainTensorVolumetric(obj,phi)
            mI = obj.materialInterpolator;
            [mu,~] = mI.computeConstitutiveTensorParams(phi);
            kappa  = ConstantFunction.create(0,obj.mesh);
            V = obj.createMaterial(mu,kappa);
        end

        function D = obtainTensorDeviatoric(obj,phi)
            mI = obj.materialInterpolator;
            [~,kappa] = mI.computeConstitutiveTensorParams(phi);
            mu        = ConstantFunction.create(0,obj.mesh);
            D = obj.createMaterial(mu,kappa);
        end

        function k = getBulk(obj,phi)
            mI = obj.materialInterpolator;
            [mu,kappa] = mI.computeConstitutiveTensorDerivativeParams(phi);
            k = obj.createMaterial(mu,kappa);
        end

        function dk = get

        function ddC = obtainTensorSecondDerivative(obj,phi)
            mI = obj.materialInterpolator;
            [mu,kappa] = mI.computeConstitutiveTensorSecondDerivativeParams(phi);
            ddC = obj.createMaterial(mu,kappa);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.materialInterpolator  = MaterialInterpolator.create(cParams.interp);
        end

    end

    methods (Access = private)

        function mat = createMaterial(obj,mu,kappa)
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = obj.mesh.ndim;
            s.shear   = mu;
            s.bulk    = kappa;
            mat = Material.create(s);
        end     

    end

    end

     %% Energy split mode
    methods (Access = public)

        function mat = getBulkMaterial(obj,u,phi,interpType)
            mu    = ConstantFunction.create(0,obj.mesh);
            kappa = obj.getBulkFun(u,phi,interpType);
            s.shear = mu;
            s.bulk  = kappa;
            s.ndim  = obj.mesh.ndim;
            mat = Isotropic2dElasticMaterial(s);
        end

        function mat = getShearMaterial(obj,phi,interpType)
            mu    = obj.getShearFun(phi,interpType);
            kappa = ConstantFunction.create(0,obj.mesh);
            s.shear = mu;
            s.bulk  = kappa;
            s.ndim  = obj.mesh.ndim;
            mat = Isotropic2dElasticMaterial(s);
        end
        
        function kFun = getBulkFun(obj,u,phi,interpType)
            fun = obj.selectDegradationFun(interpType);
            g = obj.computeSplitDegradationFun(fun,u,phi);

            E  = obj.baseMaterial.young.constant;
            nu = obj.baseMaterial.poisson.constant;
            N = obj.mesh.ndim;
            kV = E./(N*(1-(N-1)*nu));
            k = ConstantFunction.create(kV,obj.mesh);

            kFun = g.*k;
        end

        function muFun = getShearFun(obj,phi,interpType)
            fun = obj.selectDegradationFun(interpType);
            g = obj.computeDegradationFun(fun,phi);

            E  = obj.baseMaterial.young.constant;
            nu = obj.baseMaterial.poisson.constant;
            muV = E./(2*(1+nu));
            mu = ConstantFunction.create(muV,obj.mesh);

            muFun = g.*mu;
        end

    end

    methods (Access = private)

        function fun = selectDegradationFun(obj,interpType)
            switch interpType
                case 'Interpolated'
                    fun = obj.degradation.fun;
                case 'Jacobian'
                    fun = obj.degradation.dfun;
                case 'Hessian'
                    fun = obj.degradation.ddfun;
            end
        end

        function g = computeDegradationFun(obj,fun,phi)
            s.operation = @(xV) fun.evaluate(phi.evaluate(xV));
            s.ndimf = 1;
            s.mesh = obj.mesh;
            g = DomainFunction(s);
        end

        function g = computeSplitDegradationFun(obj,fun,u,phi)
            g0 = obj.computeDegradationFun(fun,phi);
            trcSign = Heaviside(trace(AntiVoigt(SymGrad(u))));
            g = g0.*trcSign + (1-trcSign);
        end
       