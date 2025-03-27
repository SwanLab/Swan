classdef MaterialPhaseFieldAnalytical < Material

    properties (Access = private)
        materialInterpolator
        mesh
    end

    methods (Access = public)

        function obj = MaterialPhaseFieldAnalytical(cParams)
            obj.init(cParams)
        end

        function C = obtainTensor(obj,phi)
            mI = obj.materialInterpolator;
            [mu,kappa] = mI.computeConstitutiveTensorParams(phi);
            C = obj.createMaterial(mu,kappa);
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

        function dC = obtainTensorDerivative(obj,phi)
            mI = obj.materialInterpolator;
            [mu,kappa] = mI.computeConstitutiveTensorDerivativeParams(phi);
            dC = obj.createMaterial(mu,kappa);
        end

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