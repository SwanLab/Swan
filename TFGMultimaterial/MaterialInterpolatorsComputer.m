classdef MaterialInterpolatorsComputer < handle

    properties (Access = public)
        materialInterpolator
    end

    properties (Access = private)
        mat
    end

    methods (Access = public)

        function obj = MaterialInterpolatorsComputer(cParams)
            obj.init(cParams);
            obj.createMaterialInterpolator1();
            obj.createMaterialInterpolator2();
            obj.createMaterialInterpolator3();
            obj.createMaterialInterpolator4();
            obj.createMaterialInterpolator5();
            obj.createMaterialInterpolator6(); 
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.mat = cParams.mat;
        end

        function createMaterialInterpolator1(obj)
            % Interpolator between VOID and MATERIAL 1
            materialA = obj.mat.D;
            materialB = obj.mat.A;
            E0   = materialA.young;
            nu0  = materialA.nu;
            E1   = materialB.young;
            nu1  = materialB.nu;
            ndim = 2;

            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator{1} = m;
        end 

        function createMaterialInterpolator2(obj)
            % Interpolator between VOID and MATERIAL 2
            materialA = obj.mat.D;
            materialB = obj.mat.B;
            E0   = materialA.young;
            nu0  = materialA.nu;
            E1   = materialB.young;
            nu1  = materialB.nu;
            ndim = 2;

            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator{2} = m;
        end 

        function createMaterialInterpolator3(obj)
            % Interpolator between VOID and MATERIAL 3
            materialA = obj.mat.D;
            materialB = obj.mat.C;
            E0   = materialA.young;
            nu0  = materialA.nu;
            E1   = materialB.young;
            nu1  = materialB.nu;
            ndim = 2;

            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator{3} = m;
        end

        function createMaterialInterpolator4(obj)
            % Interpolator between MATERIAL 1 and MATERIAL 2
            materialA = obj.mat.A;
            materialB = obj.mat.B;
            E0   = materialA.young;
            nu0  = materialA.nu;
            E1   = materialB.young;
            nu1  = materialB.nu;
            ndim = 2;

            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator{4} = m;
        end

        function createMaterialInterpolator5(obj)
            % Interpolator between MATERIAL 1 and MATERIAL 3
            materialA = obj.mat.A;
            materialB = obj.mat.C;
            E0   = materialA.young;
            nu0  = materialA.nu;
            E1   = materialB.young;
            nu1  = materialB.nu;
            ndim = 2;

            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator{5} = m;
        end

        function createMaterialInterpolator6(obj)
            % Interpolator between MATERIAL 2 and MATERIAL 3
            materialA = obj.mat.B;
            materialB = obj.mat.C;
            E0   = materialA.young;
            nu0  = materialA.nu;
            E1   = materialB.young;
            nu1  = materialB.nu;
            ndim = 2;

            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator{2} = m;
        end 
    end
end