classdef Material_Interpolation < handle
    
    properties (SetAccess = protected, GetAccess = public)
        E1
        E0
        nu1
        nu0
        rho1
        rho0        
    end
    
    properties (Access = protected)
        nstre
        ndim        
        nElem
        muFunc
        dmuFunc        
        kappaFunc
        dkappaFunc          
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
            cP = cParams.constitutiveProperties;  
            obj.rho1 = cP.rho_plus;
            obj.rho0 = cP.rho_minus;
            obj.E1   = cP.E_plus;
            obj.E0   = cP.E_minus;
            obj.nu1  = cP.nu_plus;
            obj.nu0  = cP.nu_minus;            
        end
        
        function computeSymbolicInterpolationFunctions(obj)
            [muS,dmuS,kS,dkS] = obj.computeSymbolicMuKappa();            
            obj.muFunc     = matlabFunction(muS);            
            obj.dmuFunc    = matlabFunction(dmuS);
            obj.kappaFunc  = matlabFunction(kS);
            obj.dkappaFunc = matlabFunction(dkS);            
        end
        
    end
    
    methods (Access = protected, Abstract)
        computeSymbolicMuKappa(obj)
    end
    
end