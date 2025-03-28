classdef PhaseFieldAmbrosioTortorelliSplit < handle
    
   properties (Access = private)
        shear
        bulk
        pExp
   end

    methods (Access = public)
        function obj = PhaseFieldAmbrosioTortorelliSplit(cParams)
            obj.init(cParams)
        end

        function k = computeBulkFunction(obj,u,phi)
            k0 = obj.bulk;
            k  = obj.interpolate(phi,k0);
            k  = obj.applySplit(u,k0,k);
        end

        function mu = computeShearFunction(obj,phi)
            mu0 = obj.shear;
            mu  = obj.interpolate(phi,mu0);
        end

        function dk = computeBulkFunctionDerivative(obj,u,phi)
            k0 = obj.bulk;
            dk = obj.derive(phi,k0);
            dk = obj.applySplit(u,k0,dk);
        end

        function dmu = computeShearFunctionDerivative(obj,phi)
            mu0 = obj.shear;
            dmu = obj.derive(phi,mu0);
        end
        
        function ddk = computeBulkSecondDerivative(obj,u,phi)
            k0  = obj.bulk;
            ddk = obj.derive2(phi,k0);
            ddk = obj.applySplit(u,k0,ddk);
        end

        function ddmu = computeShearFunctionSecondDerivative(obj,phi)
            mu0  = obj.shear;
            ddmu = obj.derive2(phi,mu0);
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            if isfield(cParams,'young') && isfield(cParams,'poisson')
                ndim = cParams.ndim;
                E    = cParams.young;
                nu   = cParams.poisson;
                obj.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E,nu,ndim);
                obj.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E,nu);
            elseif isfield(cParams,'bulk') && isfield(cParams,'shear')
                obj.shear = cParams.shear;
                obj.bulk  = cParams.bulk;
            end
            obj.pExp  = 2;
        end

        function f = interpolate(obj,phi,f0)
            p = obj.pExp;
            f = ((1-phi).^p).*f0;
        end
        
        function f = derive(obj,phi,f0)
            p = obj.pExp;
            f = -p*((1-phi).^(p-1)).*f0;
        end

        function f = derive2(obj,phi,f0)
            p = obj.pExp;
            f = p*(p-1)*((1-phi).^(p-2)).*f0;
        end

        function f = applySplit(obj,u,f0,f)
            trcSign = Heaviside(trace(AntiVoigt(SymGrad(u))));
            f = f.*trcSign + f0.*(1-trcSign);
        end
    end

end