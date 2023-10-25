classdef PhaseFieldInternalEnergyInterpolator < handle

    properties (Access = public)

    end

    properties (Access = private)
        isoMaterial
        mu0
        mu1
        kappa0
        kappa1
        mu
        dmu
        ddmu
        kappa
        dkappa
        ddkappa
    end

    properties (Access = private)
       constitutiveProperties        
    end

    methods (Access = public)

        function obj = PhaseFieldInternalEnergyInterpolator(cParams)
            obj.init(cParams)    
            obj.createIsotropicMaterial()
            obj.computeShearAndBulkModulus()
            obj.createInternalEnergyInterpolation()
        end

        function mat = computeMatProp(obj,phi)
            mat.mu = obj.mu(phi);
            mat.kappa = obj.kappa(phi);
        end

        function mat = computeDMatProp(obj,phi)
            mat.dmu = obj.dmu(phi);
            mat.dkappa = obj.dkappa(phi);
        end

        function mat = computeDDMatProp(obj,phi)
            mat.ddmu = obj.ddmu(phi);
            mat.ddkappa = obj.ddkappa(phi);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.constitutiveProperties = cParams.constitutiveProperties;
        end

        function createIsotropicMaterial(obj)
            s.pdim  = '2D';
            s.ptype = 'ELASTIC';
            obj.isoMaterial = Material.create(s);
        end        

        function computeShearAndBulkModulus(obj)
             E1 = obj.constitutiveProperties.E_plus;
             nu1 = obj.constitutiveProperties.nu_plus;
             obj.mu1 = obj.isoMaterial.computeMuFromYoungAndNu(E1,nu1);
             obj.kappa1 = obj.isoMaterial.computeKappaFromYoungAndNu(E1,nu1);

             E0 = obj.constitutiveProperties.E_minus;
             nu0 = obj.constitutiveProperties.nu_minus;
             obj.mu0 = obj.isoMaterial.computeMuFromYoungAndNu(E0,nu0);
             obj.kappa0 = obj.isoMaterial.computeKappaFromYoungAndNu(E0,nu0);
        end

        function createInternalEnergyInterpolation(obj)
            m0 = obj.mu0;
            m1 = obj.mu1;
            k0 = obj.kappa0;
            k1 = obj.kappa1;
            p = 2;

            obj.mu      = @(phi) (1-phi).^p*m1 + m0;
            obj.kappa   = @(phi) (1-phi).^p*k1 + k0;
            obj.dmu     = @(phi) -p*(1-phi).^(p-1)*m1 + m0;
            obj.dkappa  = @(phi) -p*(1-phi).^(p-1)*k1 + k0;
            obj.ddmu    = @(phi) p*(p-1)*(1-phi).^(p-2)*m1 + m0;
            obj.ddkappa = @(phi) p*(p-1)*(1-phi).^(p-2)*k1 + k0;
        end

    end

end