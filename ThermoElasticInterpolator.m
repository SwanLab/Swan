classdef ThermoElasticInterpolator < handle
    
    properties (Access = public)
        fun
        dfun
    end

    properties (Access = private)
        mesh
        thermalProblem
        elasticInt
    end

    methods (Access = public)
        
        function obj = ThermoElasticInterpolator(cParams)
            obj.init(cParams)
        end

        function [mu,kappa] = computeConsitutiveTensor(obj,rho)
            % 1. Interpolate kappa with density
            % 2. Solve thermal problem
            % 3. Create elastic interpolator

            mu    = obj.computeMuFunction(rho{1});
            kappa = obj.computeKappaFunction(rho{1});
        end

        function [dmu,dkappa] = computeConsitutiveTensorDerivative(obj,rho)
            dmu{1}    = [];
            dkappa{1} = [];
        end

    end

    methods (Access = private)
        
        function init(obj,cParams) % Composition with thermal problem + elastic interpolator
            obj.mesh      = cParams.mesh;
            obj.thermalProblem = cParams.thermalProblem;
        end

        function  createElasticInterpolator(obj,s,T)
            alpha = 4e-4;

            E0 = s.E0;
            nu0 = s.nu0;
            ndim = obj.mesh.ndim;
            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            E1 = 1*(1-alpha.*T);
            nu1 = s.nu1;
            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.interpolation  = 'SIMPALL';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.elasticInt = m;
        end

    end

end