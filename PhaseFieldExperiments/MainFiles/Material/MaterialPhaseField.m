classdef MaterialPhaseField < Material

    properties (Access = private)
        phi
        u
    end

    properties (Access = private)
        mesh
        degradation
        baseMaterial
        Gc
    end

    methods (Access = public)

        function obj = MaterialPhaseField(cParams)
            obj.init(cParams)
        end

        function obj = setDesignVariable(obj,u,phi)
            obj.u = u;
            obj.phi = phi;
        end

        function C = obtainTensor(obj)
            f    = obj.degradation.fun;
            C{1} = obj.createDegradedMaterial(f);
        end

        function dC = obtainTensorDerivative(obj)
            df    = obj.degradation.dfun;
            dC{1} = obj.createDegradedMaterial(df);
        end

        function ddC = obtainTensorSecondDerivative(obj)
            ddf    = obj.degradation.ddfun;
            ddC{1} = obj.createDegradedMaterial(ddf);
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

        function mat = createBaseMaterial(obj,cParams)
            sIso.ndim = obj.mesh.ndim;
            sIso.young = ConstantFunction.create(cParams.matInfo.E,obj.mesh);
            sIso.poisson = ConstantFunction.create(cParams.matInfo.nu,obj.mesh);
            mat = Isotropic2dElasticMaterial(sIso);
        end

        function mat = createDegradedMaterial(obj,degFun)
            df    = degFun;
            mu    = obj.baseMaterial.createShear();
            kappa = obj.baseMaterial.createBulk();
            degM  = obj.createDegradedLameParameterFunction(mu,df);
            degK  = obj.createDegradedLameParameterFunction(kappa,df);
            s.shear = degM;
            s.bulk  = degK;
            s.ndim  = obj.mesh.ndim;
            mat = Isotropic2dElasticMaterial(s);
        end     
    
        function xf = createDegradedLameParameterFunction(obj,param,f)
            s.operation = @(xV) param.evaluate(xV).*obj.evaluateDegradation(f,xV);
            s.ndimf = 1;
            xf =  DomainFunction(s);
        end      

        function fV = evaluateDegradation(obj,f,xV)
            phiV = obj.phi.evaluate(xV);
            fV = f(phiV);
        end            


    end
    
    %% Energy split mode
    methods (Access = public)
        
        function kFun = getBulkFun(obj,u,phi,interpType)
            obj.setDesignVariable(u,phi);
            [~,g0] = obj.computeDegradationFun(interpType);

            E  = obj.baseMaterial.young.constant;
            nu = obj.baseMaterial.poisson.constant;
            N = 2;
            kV = E./(N*(1-(N-1)*nu));
            k = ConstantFunction.create(kV,obj.mesh);

            kFun = g0.*k';
        end

        function muFun = getShearFun(obj,u,phi,interpType)
            obj.setDesignVariable(u,phi);
            [~,g0] = obj.computeDegradationFun(interpType);

            E  = obj.baseMaterial.young.constant;
            nu = obj.baseMaterial.poisson.constant;
            muV = E./(2*(1+nu));
            mu = ConstantFunction.create(muV,obj.mesh);

            muFun = g0.*mu;
        end

        function mat = getBulkMaterial(obj,u,phi)
            obj.setDesignVariable(u,phi);
            df    = obj.degradation.fun;
            mu    = ConstantFunction.create(0,obj.mesh);
            kappa = obj.baseMaterial.createBulk();
            degM  = obj.createDegradedLameParameterFunction(mu,df);
            degK  = obj.createDegradedLameParameterFunction(kappa,df);
            s.shear = degM;
            s.bulk  = degK;
            s.ndim  = obj.mesh.ndim;
            mat = Isotropic2dElasticMaterial(s);
        end

        function mat = getShearMaterial(obj,u,phi)
            obj.setDesignVariable(u,phi);
            df    = obj.degradation.fun;
            mu    = obj.baseMaterial.createShear();
            kappa = ConstantFunction.create(0,obj.mesh);
            degM  = obj.createDegradedLameParameterFunction(mu,df);
            degK  = obj.createDegradedLameParameterFunction(kappa,df);
            s.shear = degM;
            s.bulk  = degK;
            s.ndim  = obj.mesh.ndim;
            mat = Isotropic2dElasticMaterial(s);
        end
    end

    methods (Access = private)

        function [g, g0] = computeDegradationFun(obj,interpType)
            switch interpType
                case 'Interpolated'
                    fun = obj.degradation.fun;
                case 'Jacobian'
                    fun = obj.degradation.dfun;
                case 'Hessian'
                    fun = obj.degradation.ddfun;
            end
            s.operation = @(xV) fun(obj.phi.evaluate(xV));
            s.ndimf = obj.phi.ndimf;
            g0 = DomainFunction(s);

            trcSign = Heaviside(trace(SymGrad(obj.u)));
            g = g0.*trcSign + (1-trcSign);
        end

    end



end