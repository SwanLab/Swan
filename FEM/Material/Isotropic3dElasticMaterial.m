classdef Isotropic3dElasticMaterial < IsotropicElasticMaterial
    
    methods (Access = public)
        
        function obj = Isotropic3dElasticMaterial(cParams)
            obj.init(cParams);
            obj.nstre = 6;
        end
        
    end
    
    methods (Access = protected)

        function computeC(obj)   
            m = obj.mu;
            l = obj.computeLambdaFromMuAndKappa(obj.mu,obj.kappa);
            C = zeros(obj.nstre,obj.nstre,obj.nElem,obj.nGaus);            
            C(1,1,:,:) = 2*m + l;
            C(2,2,:,:) = 2*m + l;
            C(3,3,:,:) = 2*m + l;            
            C(1,2,:,:) = l;
            C(2,1,:,:) = l;
            C(1,3,:,:) = l;
            C(3,1,:,:) = l;
            C(3,2,:,:) = l;
            C(2,3,:,:) = l;
            C(4,4,:,:) = m;
            C(5,5,:,:) = m;
            C(6,6,:,:) = m;
            obj.C = C;
        end
        
    end
    
    methods (Access = public, Static)
        
        function E = computeYoungFromMuAndKappa(mu,kappa)
            E = 9*kappa*mu/(3*kappa + mu);
        end
        
        function nu = computeNuFromMuAndKappa(mu,kappa)
            nu = (3*kappa-2*mu)/(6*kappa+2*mu);
        end        
        
        function k = computeKappaFromYoungAndNu(E,nu)
            k = E/(3*(1-2*nu));
        end
        
        function lambda = computeLambdaFromMuAndKappa(mu,kappa)
            lambda = kappa - 2/3*mu;
        end
        
        function kappa = computeKappaFromMuAndLambda(mu,lambda)
            kappa = lambda + 2/3*mu;
        end
        
    end
end

