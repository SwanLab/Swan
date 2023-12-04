classdef PhaseFieldInternalEnergyInterpolator < handle

    properties (Access = public)

    end

    properties (Access = private)
        isoMaterial
        g
        d1g
        d2g
        mu
        kappa
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

        function g = computeDegradation(obj,phi,deriv)
            if deriv == 0
                g = obj.g(phi);
            elseif deriv == 1
                g = obj.d1g(phi);
            elseif deriv == 2
                g = obj.d2g(phi);
            else
                error(['Degree ',num2str(deriv),' derivative not implemented'])
            end
        end

        function mat = computeMatProp(obj,phi,deriv)
            m = obj.mu;
            k = obj.kappa;
            if deriv == 0
                mat.mu    = obj.g(phi)*m;
                mat.kappa = obj.g(phi)*k;
            elseif deriv == 1
                mat.mu    = obj.d1g(phi)*m;
                mat.kappa = obj.d1g(phi)*k;
            elseif deriv == 2
                mat.mu    = obj.d2g(phi)*m;
                mat.kappa = obj.d2g(phi)*k;
            else
                error(['Degree ',num2str(deriv),' derivative not implemented'])
            end
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
             E = obj.constitutiveProperties.E;
             nu = obj.constitutiveProperties.nu;
             obj.mu = obj.isoMaterial.computeMuFromYoungAndNu(E,nu);
             obj.kappa = obj.isoMaterial.computeKappaFromYoungAndNu(E,nu);
        end

        function createInternalEnergyInterpolation(obj)
            p = 2;

            obj.g      = @(phi) ((1-phi).^p);
            obj.d1g     = @(phi) -p*((1-phi).^(p-1));
            obj.d2g = @(phi) p*(p-1)*((1-phi).^(p-2));
        end

    end

end