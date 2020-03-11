classdef Isotropic3dElasticMaterial < IsotropicElasticMaterial
    
    methods (Access = public)
        
        function obj = Isotropic3dElasticMaterial(cParams)
            obj.init(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function computeLambda(obj)
            obj.lambda = obj.kappa - 2/3*obj.mu;
        end
        
        function obj = computeC(obj)
            m = obj.mu;
            l = obj.lambda;
            obj.C(1,1,:) = 2*m + l;
            obj.C(2,2,:) = 2*m + l;
            obj.C(3,3,:) = 2*m + l;            
            obj.C(1,2,:) = l;
            obj.C(2,1,:) = l;
            obj.C(1,3,:) = l;
            obj.C(3,1,:) = l;
            obj.C(3,2,:) = l;
            obj.C(2,3,:) = l;
            obj.C(4,4,:) = m;
            obj.C(5,5,:) = m;
            obj.C(6,6,:) = m;
        end
        
    end
    
%     methods (Access = private)
%         
%         function computeYoungModulus(obj)
%             k = obj.kappa;
%             m = obj.mu;
%             E = 9*k.*m./(3*k + m);
%             obj.E = E;
%         end
%         
%         function computePoissonRatio(obj)
%             k = obj.kappa;
%             m = obj.mu;
%             nu = (3*k - 2*m)./(2*(3*k + m));
%             obj.nu = nu;
%         end
%         
%     end
    
end

