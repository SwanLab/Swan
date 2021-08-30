classdef Isotropic2dElasticMaterial < IsotropicElasticMaterial

    methods (Access = public)
        
        function obj = Isotropic2dElasticMaterial(cParams)
            obj.init(cParams);    
            obj.nstre = 3;            
        end
        
    end    
    
    methods (Access = protected)
        
        function C = computeC(obj)
            m = obj.mu;
            l = obj.computeLambdaFromMuAndKappa(obj.mu,obj.kappa);
            C = zeros(obj.nstre,obj.nstre,obj.nElem,obj.nGaus);                        
            C(1,1,:,:)= 2*m+l;
            C(1,2,:,:)= l;
            C(2,1,:,:)= l;
            C(2,2,:,:)= 2*m+l;
            C(3,3,:,:)= m;
            obj.C = C;
        end
        
    end
    
    methods (Access = public, Static)
        
        function E = computeYoungFromMuAndKappa(mu,kappa)
            E = 4*kappa*mu/(kappa + mu);
        end
        
        function nu = computeNuFromMuAndKappa(mu,kappa)
            nu = (kappa - mu)/(kappa +mu);
        end        
        
        function k = computeKappaFromYoungAndNu(E,nu)
            k = E/(2*(1-nu));
        end
        
        function lambda = computeLambdaFromMuAndKappa(mu,kappa)
            lambda = kappa - mu;
        end
        
        function kappa = computeKappaFromMuAndLambda(mu,lambda)
            kappa = lambda + mu;
        end
        
    end

end

