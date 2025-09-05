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
            if isa(mI,'SimpAllExplicitInterpolator')
                rho = (1-phi.fun);
                [mu,kappa] = mI.computeConstitutiveTensor(rho);
            else
                [mu,kappa] = mI.computeConstitutiveTensor(phi);
            end
            C = obj.createMaterial(mu,kappa);
        end

        function dC = obtainTensorDerivative(obj,phi)
            mI = obj.materialInterpolator;
            if isa(mI,'SimpAllExplicitInterpolator')
                rho = (1-phi.fun);
                [dmu,dkappa] = mI.computeConstitutiveTensorDerivative(rho);
                dmu = -dmu; dkappa = -dkappa;
            else
                [dmu,dkappa] = mI.computeConstitutiveTensorDerivative(phi);
            end
            dC = obj.createMaterial(dmu,dkappa);
        end

        function ddC = obtainTensorSecondDerivative(obj,phi)
            mI = obj.materialInterpolator;
            if isa(mI,'SimpAllExplicitInterpolator')
                rho = (1-phi.fun);
                [ddmu,ddkappa] = mI.computeConstitutiveTensorDerivative(rho); %Just to compute something
            else
                [ddmu,ddkappa] = mI.computeConstitutiveTensorSecondDerivative(phi);
            end
            ddC = obj.createMaterial(ddmu,ddkappa);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            if isfield(cParams.interp,'subType')
                if cParams.interp.subType == "ATSplit"
                    fprintf('ATSplit not supported by this material, changing to AT \n')
                    cParams.interp.subType = 'AT';
                end
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