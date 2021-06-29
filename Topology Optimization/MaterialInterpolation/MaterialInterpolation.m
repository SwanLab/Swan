classdef MaterialInterpolation < handle
    
   properties (Access = protected)
        nstre
        ndim    
        pdim
        nElem
        muFunc
        dmuFunc        
        kappaFunc
        dkappaFunc      
            
        matProp
        isoMaterial
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = MaterialInterpolationFactory;
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public)
                
        function mp = computeMatProp(obj,rho)
            mu      = obj.muFunc(rho);
            kappa   = obj.kappaFunc(rho);
            dmu     = obj.dmuFunc(rho);
            dkappa  = obj.dkappaFunc(rho);            
            mp.mu     = mu;
            mp.kappa  = kappa;
            mp.dmu    = dmu;
            mp.dkappa = dkappa;
        end                           
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.nElem = cParams.nElem;  
            obj.ndim  = cParams.ndim;
            obj.pdim  = cParams.dim;
            obj.saveYoungAndPoisson(cParams);
            obj.createIsotropicMaterial();            
            obj.computeMuAndKappaIn0();
            obj.computeMuAndKappaIn1();
        end
        
        function computeMuAndKappaIn0(obj)
            E0  = obj.matProp.E0;
            nu0 = obj.matProp.nu0;
            obj.matProp.mu0    = obj.computeMu(E0,nu0);
            obj.matProp.kappa0 = obj.computeKappa(E0,nu0);
           
        end
        
        function computeMuAndKappaIn1(obj)
            E1  = obj.matProp.E1;
            nu1 = obj.matProp.nu1;            
            obj.matProp.mu1    = obj.computeMu(E1,nu1);            
            obj.matProp.kappa1 = obj.computeKappa(E1,nu1);              
        end
        
        function saveYoungAndPoisson(obj,cParams)
            cP = cParams.constitutiveProperties;  
            obj.matProp.rho1 = cP.rho_plus;
            obj.matProp.rho0 = cP.rho_minus;
            obj.matProp.E1   = cP.E_plus;
            obj.matProp.E0   = cP.E_minus;
            obj.matProp.nu1  = cP.nu_plus;
            obj.matProp.nu0  = cP.nu_minus;            
        end
        
        function createIsotropicMaterial(obj)
            s.pdim  = obj.pdim;            
            s.ptype = 'ELASTIC';
            obj.isoMaterial = Material.create(s);                      
        end
        
        function computeSymbolicInterpolationFunctions(obj)
            [muS,dmuS,kS,dkS] = obj.computeSymbolicMuKappa();            
            obj.muFunc     = matlabFunction(muS);            
            obj.dmuFunc    = matlabFunction(dmuS);
            obj.kappaFunc  = matlabFunction(kS);
            obj.dkappaFunc = matlabFunction(dkS);            
        end
        
        function [muS,dmuS,kS,dkS] = computeSymbolicMuKappa(obj)
            [muS,dmuS] = obj.computeMuSymbolicFunctionAndDerivative();
            [kS,dkS]   = obj.computeKappaSymbolicFunctionAndDerivative();   
        end    
        
        function mu = computeMu(obj,E,nu)
            mat = obj.isoMaterial;
            mu = mat.computeMuFromYoungAndNu(E,nu);    
        end
        
        function k = computeKappa(obj,E,nu)
            mat = obj.isoMaterial;            
            k = mat.computeKappaFromYoungAndNu(E,nu);    
        end       
        
    end
    
    methods (Access = protected, Abstract)
        computeMuSymbolicFunctionAndDerivative(obj)
        computeKappaSymbolicFunctionAndDerivative(obj)        
    end
    
end