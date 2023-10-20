classdef PhaseFieldInterpolator < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        isoMaterial
        mu0
        mu1
        kappa1
        kappa0
    end
    
    properties (Access = private)
       constitutiveProperties        
    end
    
    methods (Access = public)
        
        function obj = PhaseFieldInterpolator(cParams)
            obj.init(cParams)    
            obj.createIsotropicMaterial()
            obj.computeShearAndBulkModulus()
        end
        
        function mat = computeMatProp(obj,phi)
            % Constant distribution
            mat.mu = obj.mu1;
            mat.kappa = obj.kappa1;

            % Linear distribution of g(phi)
            % mat.mu = obj.mu1*phi + obj.mu0*(1-phi);
            % mat.kappa = obj.kappa1*phi + obj.kappa0*(1-phi);
            
            % Quadratic distribution of g(phi)
            % mat.mu = (obj.mu1-obj.mu0)*phi.^2+ obj.mu0;
            % mat.kappa = (obj.kappa1-obj.kappa0)*phi.^2+ obj.kappa0;
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
        
    end
    
end