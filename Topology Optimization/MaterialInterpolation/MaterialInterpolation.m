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


        % -------------- SIMP modal ------------------------------
        muLFunc
        muHFunc
        dmuLFunc
        dmuHFunc
        kappaLFunc
        kappaHFunc
        dkappaLFunc
        dkappaHFunc
        % --------------------------------------------------------
   end

   methods (Access = public, Static)
        
        function obj = create(cParams)
            f = MaterialInterpolationFactory;
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public)
        
        function mp = computeMatProp(obj,rho)
%             mu      = obj.muFunc(rho);
%             kappa   = obj.kappaFunc(rho);
%             dmu     = obj.dmuFunc(rho);
%             dkappa  = obj.dkappaFunc(rho);

            %-------------- SIMP_modal ----------------------
            mu = zeros(length(rho),1);
            kappa = zeros(length(rho),1);
            dmu = zeros(length(rho),1);
            dkappa = zeros(length(rho),1);
            for iElem = 1:length(rho)
                if rho(iElem) <= 0.1
                    mu(iElem,1)      = obj.muLFunc(rho(iElem));
                    kappa(iElem,1)   = obj.kappaLFunc(rho(iElem));
                    dmu(iElem,1)     = obj.dmuLFunc(); %7.4925e-3; % 3.74625e-4;% 1.873125e-3 ;% obj.dmuLFunc(rho(iElem));
                    dkappa(iElem,1)  = obj.dkappaLFunc(); % 7.4925e-4; % 3.74625e-3 ;%obj.dkappaLFunc(rho(iElem));
                else
                    mu(iElem,1)      = obj.muHFunc(rho(iElem));
                    kappa(iElem,1)   = obj.kappaHFunc(rho(iElem));
                    dmu(iElem,1)     = obj.dmuHFunc(rho(iElem));
                    dkappa(iElem,1)  = obj.dkappaHFunc(rho(iElem));
                end
            end
%           %------------------------------------

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

        % -------------- SIMP modal ------------------------------
        function computeSymbolicInterpolationFunctions(obj)
            [muLS,muHS,dmuLS,dmuHS,kLS,kHS,dkLS,dkHS] = obj.computeSymbolicMuKappa();
            obj.muLFunc     = matlabFunction(muLS);
            obj.muHFunc     = matlabFunction(muHS);
            obj.dmuLFunc    = matlabFunction(dmuLS);
            obj.dmuHFunc    = matlabFunction(dmuHS);
            obj.kappaLFunc  = matlabFunction(kLS);
            obj.kappaHFunc  = matlabFunction(kHS);
            obj.dkappaLFunc = matlabFunction(dkLS);
            obj.dkappaHFunc = matlabFunction(dkHS);
        end
         % -----------------------------------------------------

        % -------------- SIMP modal ------------------------------
        function [muLS,muHS,dmuLS,dmuHS,kLS,kHS,dkLS,dkHS] = computeSymbolicMuKappa(obj)
            [muLS,muHS,dmuLS,dmuHS] = obj.computeMuSymbolicFunctionAndDerivative();
            [kLS,kHS,dkLS,dkHS]   = obj.computeKappaSymbolicFunctionAndDerivative;
        end
         % -----------------------------------------------------
%         function computeSymbolicInterpolationFunctions(obj)
%             [muS,dmuS,kS,dkS] = obj.computeSymbolicMuKappa();
%             obj.muFunc     = matlabFunction(muS);
%             obj.dmuFunc    = matlabFunction(dmuS);
%             obj.kappaFunc  = matlabFunction(kS);
%             obj.dkappaFunc = matlabFunction(dkS);
%         end
        
%         function [muS,dmuS,kS,dkS] = computeSymbolicMuKappa(obj)
%             [muS,dmuS] = obj.computeMuSymbolicFunctionAndDerivative();
%             [kS,dkS]   = obj.computeKappaSymbolicFunctionAndDerivative;
%         end
        
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