classdef MaterialPropertiesComputer < handle

    properties (Access = public)
        matA
        matB
        matC
        matD
        ndim
    end

    methods (Access = public)

        function obj = MaterialPropertiesComputer()
            obj.init();
            obj.createMaterialA();
            obj.createMaterialB();
            obj.createMaterialC();
            obj.createMaterialD();
        end
    end

    methods (Access = private)

        function init(obj)
            obj.matA.young = 200E9;
            obj.matA.nu = 0.25;
            
            obj.matB.young = 100E9;
            obj.matB.nu = 0.25;
           
            obj.matC.young = 50E9;
            obj.matC.nu = 0.25;
            
            obj.matD.young = 200E9*1/1000;
            obj.matD.nu = 0.25;

            obj.ndim = 2;
        end

        function createMaterialA(obj)
            E0 = obj.matA.young;
            nu0 = obj.matA.nu;
            dim = obj.ndim;

            %obj.matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            la = nu0.*E0./((1+nu0).*(1-2.*nu0)); % plain strain
            obj.matA.shear = E0./(2.*(1+nu0));
            mu0 = obj.matA.shear;
            obj.matA.lambda = 2.*mu0.*la./(la+2.*mu0); % plane stress
            %obj.matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,dim);
            %obj.matA.lambda = IsotropicElasticMaterial.computeLambdaFromShearAndBulk(E0,nu0,dim);
        end

        function createMaterialB(obj)
            E1 = obj.matB.young;
            nu1 = obj.matB.nu;
            dim = obj.ndim;

            %obj.matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            la = nu1.*E1./((1+nu1).*(1-2.*nu1)); 
            obj.matB.shear = E1./(2.*(1+nu1));
            mu1 = obj.matB.shear;
            obj.matB.lambda = 2.*mu1.*la./(la+2.*mu1);
            obj.matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,dim);
            %obj.matB.lambda = IsotropicElasticMaterial.computeLambdaFromShearAndBulk(E1,nu1,dim);
        end

        function createMaterialC(obj)
            E2 = obj.matC.young;
            nu2 = obj.matC.nu;
            dim = obj.ndim;

            %obj.matC.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E2,nu2);
            la = nu2.*E2./((1+nu2).*(1-2.*nu2)); 
            obj.matC.shear = E2./(2.*(1+nu2));
            mu2 = obj.matC.shear;
            obj.matC.lambda = 2.*mu2.*la./(la+2.*mu2);
            obj.matC.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E2,nu2,dim);
            %obj.matC.lambda = IsotropicElasticMaterial.computeLambdaFromShearAndBulk(E2,nu2,dim);
        end

        function createMaterialD(obj)
            E3 = obj.matD.young;
            nu3 = obj.matD.nu;
            dim = obj.ndim;

            %obj.matD.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E3,nu3);
            la = nu3.*E3./((1+nu3).*(1-2.*nu3)); 
            obj.matD.shear = E3./(2.*(1+nu3));
            mu3 = obj.matD.shear;
            obj.matD.lambda = 2.*mu3.*la./(la+2.*mu3);
            obj.matD.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E3,nu3,dim);
            %obj.matD.lambda = IsotropicElasticMaterial.computeLambdaFromShearAndBulk(E3,nu3,dim);
        end
    end
end