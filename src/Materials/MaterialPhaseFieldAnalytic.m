classdef MaterialPhaseFieldAnalytic < Material

    properties (Access = private)
        materialInterpolator
        mesh
    end

    methods (Access = public)

        function obj = MaterialPhaseFieldAnalytic(cParams)
            obj.init(cParams)
        end

        function C = obtainTensor(obj,phi)
            mI = obj.materialInterpolator;
            [mu,kappa] = mI.computeConstitutiveTensorParams(phi);
            C = obj.createMaterial(mu,kappa);
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
            if cParams.interp.degFunType == "ATSplit"
                fprintf('ATSplit not supported by this material, changing to AT \n')
                cParams.interp.degFunType = 'AT';
            end
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