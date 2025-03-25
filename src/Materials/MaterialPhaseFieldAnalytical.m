classdef MaterialPhaseFieldAnalytical < Material

    properties (Access = private)
        mesh
        degradation
        baseMaterial
        Gc
    end

    methods (Access = public)

        function obj = MaterialPhaseFieldAnalytical(cParams)
            obj.init(cParams)
        end

        function C = obtainTensor(obj,phi)
            f    = obj.degradation.fun;
            degFun = obj.computeDegradationFun(f,phi);
            C = obj.createDegradedMaterial(degFun);
        end

        function dC = obtainTensorDerivative(obj,phi)
            df    = obj.degradation.dfun;
            degFun = obj.computeDegradationFun(df,phi);
            dC = obj.createDegradedMaterial(degFun);
        end

        function ddC = obtainTensorSecondDerivative(obj,phi)
            ddf    = obj.degradation.ddfun;
            degFun = obj.computeDegradationFun(ddf,phi);
            ddC = obj.createDegradedMaterial(degFun);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.degradation  = MaterialInterpolator.create(cParams);
            obj.baseMaterial = obj.createBaseMaterial(cParams);
            obj.Gc           = cParams.Gc;
        end

    end

    methods (Access = private)

        function mat = createBaseMaterial(~,cParams)
            sIso.mesh = cParams.mesh;
            sIso.ndim = cParams.mesh.ndim;
            sIso.young = ConstantFunction.create(cParams.matInfo.E,cParams.mesh);
            sIso.poisson = ConstantFunction.create(cParams.matInfo.nu,cParams.mesh);
            mat = Isotropic2dElasticMaterial(sIso);
        end

        function mat = createDegradedMaterial(obj,fun)
            E = obj.baseMaterial.young;
            nu = obj.baseMaterial.poisson;
            ndim = obj.mesh.ndim;

            mu    = obj.baseMaterial.computeMuFromYoungAndPoisson(E,nu);
            kappa = obj.baseMaterial.computeKappaFromYoungAndPoisson(E,nu,ndim);
            degM  = fun.*mu + 1e-10; 
            degK  = fun.*kappa + 1e-10;
            s.shear = degM;
            s.bulk  = degK;
            s.ndim  = obj.mesh.ndim;
            s.mesh  = obj.mesh;
            mat = Isotropic2dElasticMaterial(s);
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

    end

end