classdef MaterialPhaseField < Material
    
    properties (Access = public)
        ndimf
    end

    properties (Access = private)
        phi
        u
        tensorType
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

        function C = obtainTensor(obj)
            f    = obj.degradation.fun;
            C{1} = obj.createDegradedMaterial(f);
        end

        function dC = obtainTensorDerivative(obj)
            df    = obj.degradation.dfun;
            dC{1} = obj.createDegradedConstitutiveTensor(df);
        end

        function ddC = obtainTensorSecondDerivative(obj)
            ddf    = obj.degradation.ddfun;
            ddC{1} = obj.createDegradedConstitutiveTensor(ddf);
        end



        % function kFun = getBulkFun(obj,u,phi,interpType)
        %     obj.setDesignVariable(u,phi,interpType,'');
        %     [g,~] = obj.computeDegradationFun();
        % 
        %     E  = obj.baseMaterial.young.fValues;
        %     nu = obj.baseMaterial.poisson.fValues;
        %     kV = obj.computeKappaFromYoungAndPoisson(E,nu,obj.ndim);
        % 
        %     s.mesh = obj.mesh;
        %     s.order = obj.young.order;
        %     s.fValues = kV;
        %     k = LagrangianFunction(s);
        % 
        %     kFun = times(g,k);
        % end
        % 
        % function muFun = getShearFun(obj,u,phi,interpType)
        %     obj.setDesignVariable(u,phi,interpType,'');
        %     [~,g0] = obj.computeDegradationFun();
        % 
        %     E  = obj.baseMaterial.young.fValues;
        %     nu = obj.poisson.fValues;
        %     muV = obj.computeMuFromYoungAndPoisson(E,nu);
        % 
        %     s.mesh = obj.mesh;
        %     s.order = 'P1';
        %     s.fValues = muV;
        %     mu = LagrangianFunction(s);
        % 
        %     muFun = g0.*mu;
        % end

        function obj = setDesignVariable(obj,u,phi,tensorType)
            obj.u = u;
            obj.phi = phi;
            obj.tensorType = tensorType;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.degradation  = cParams.degradation;
            obj.baseMaterial = cParams.baseMaterial;
            obj.Gc           = cParams.Gc;
            obj.ndimf        = 9;
        end

    end

    methods (Access = private)

        function C = createDegradedConstitutiveTensor(obj,f)
            s.operation = @(xV) obj.evaluateConstiutiveTensor(xV,f);
            s.ndimf = 9;
            C =  DomainFunction(s);
        end

        function C = evaluateConstiutiveTensor(obj,xV,fun)
            dMat = obj.createDegradedMaterial(fun);
            C = dMat.evaluate(xV);
        end

        function m = createDegradedMaterial(obj,degFun)
            df    = degFun;
            mu    = obj.baseMaterial.createShear();
            kappa = obj.baseMaterial.createBulk();
            degM  = obj.createDegradedLameParameterFunction(mu,df);
            degK  = obj.createDegradedLameParameterFunction(kappa,df);
            s.shear = degM;
            s.bulk  = degK;
            s.ndim  = obj.mesh.ndim;
            m = Isotropic2dElasticMaterial(s);
        end

        % function [mu,l] = computeMatTensorParams(obj,xV,fun)
        %     [g, g0] = obj.computeDegradationFun(fun);
        %     [mu, k] = obj.computeMuAndKappa(g,g0,xV);
        %     l = obj.computeLambdaFromShearAndBulk(mu,k,obj.ndim);
        % end
        % 
        % function [g, g0] = computeDegradationFun(obj,fun)
        %     s.operation = @(xV) fun(obj.phi.evaluate(xV));
        %     s.ndimf = obj.phi.ndimf;
        %     g0 = DomainFunction(s);
        % 
        %     trcSign = Heaviside(trace(SymGrad(obj.u)));
        %     g = g0.*trcSign + (1-trcSign);
        % end
        % 
        % function [mu, k] = computeMuAndKappa(obj,g,g0,xV)
        %     [muV,kV] = obj.baseMaterial.computeShearAndBulk(xV);
        %     if obj.tensorType == "Deviatoric"
        %         kV = 0.*kV;
        %     elseif obj.tensorType == "Volumetric"
        %         muV = 0.*muV; 
        %     end
        %     mu = g0.evaluate(xV).*muV;
        %     k  = g.evaluate(xV).*kV;
        % end        
    
        function xf = createDegradedLameParameterFunction(obj,x,f)
            s.operation = @(xV) x.evaluate(xV).*obj.evaluateDegradation(f,xV);
            s.ndimf = 1;
            xf =  DomainFunction(s);
        end      

        function fV = evaluateDegradation(obj,f,xV)
            phiV = obj.phi.evaluate(xV);
            fV = f(phiV);
        end            


    end



end