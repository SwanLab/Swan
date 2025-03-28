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
            [~,kappa] = mI.computeConstitutiveTensorParams(phi);
            mu        = ConstantFunction.create(0,obj.mesh);
            V = obj.createMaterial(mu,kappa);
        end

        function D = obtainTensorDeviatoric(obj,phi,u)
            mI = obj.materialInterpolator;
            [mu,~] = mI.computeConstitutiveTensorParams(phi,u);
            kappa  = ConstantFunction.create(0,obj.mesh);
            D = obj.createMaterial(mu,kappa);
        end

        function kappa = obtainBulkFunction(obj,phi,u)
            mI    = obj.materialInterpolator;
            kappa = mI.computeBulkFunction(phi,u);
        end

        function dkappa = obtainBulkDerivative(obj,phi,u)
            mI     = obj.materialInterpolator;
            dkappa = mI.computeBulkFunctionDerivative(phi,u);
        end

        function ddkappa = obtainBulkSecondDerivative(obj,phi,u)
            mI      = obj.materialInterpolator;
            ddkappa = mI.computeBulkSecondDerivative(phi,u);
        end

        function mu = obtainShearFunction(obj,phi)
            mI = obj.materialInterpolator;
            mu = mI.computeShearFunction(phi);
        end

        function dmu = obtainShearDerivative(obj,phi)
            mI  = obj.materialInterpolator;
            dmu = mI.computeShearFunctionDerivative(phi); 
        end

        function ddmu = obtainShearSecondDerivative(obj,phi)
            mI   = obj.materialInterpolator;
            ddmu = mI.computeShearFunctionSecondDerivative(phi);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            cParams.interp.degradationType = 'ATSplit';
            obj.materialInterpolator  = MaterialInterpolator.create(cParams.interp);
        end

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
       