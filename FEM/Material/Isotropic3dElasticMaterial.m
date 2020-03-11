classdef Isotropic3dElasticMaterial < IsotropicElasticMaterial
    
    methods (Access = public)
        
        function obj = Isotropic3dElasticMaterial(cParams)
            obj.init(cParams);
            obj.ndim = 3;
        end
        
    end
    
    methods (Access = protected)
        
        function computeLambda(obj)
            d = obj.ndim;
            obj.lambda = obj.kappa - 2/d*obj.mu;
        end
        
        function computeC(obj,mu,lambda)            
            m = mu;
            l = lambda;
            C = zeros(obj.nstre,obj.nstre,obj.nElem);            
            C(1,1,:) = 2*m + l;
            C(2,2,:) = 2*m + l;
            C(3,3,:) = 2*m + l;            
            C(1,2,:) = l;
            C(2,1,:) = l;
            C(1,3,:) = l;
            C(3,1,:) = l;
            C(3,2,:) = l;
            C(2,3,:) = l;
            C(4,4,:) = m;
            C(5,5,:) = m;
            C(6,6,:) = m;
            obj.C = C;
        end
        
    end
end

