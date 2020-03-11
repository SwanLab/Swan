classdef Isotropic2dElasticMaterial < IsotropicElasticMaterial

    methods (Access = public)
        
        function obj = Isotropic2dElasticMaterial(cParams)
            obj.init(cParams);          
        end
        
    end    
    
    methods (Access = protected)
        
        function computeLambda(obj)
            obj.lambda = obj.kappa-obj.mu;                        
        end
        
        function obj = computeC(obj)
%            obj.computeYoungModulus();
%            obj.computePoissonRatio();
            m = obj.mu;
            l = obj.lambda;
            obj.C(1,1,:)= 2*m+l;
            obj.C(1,2,:)= l;
            obj.C(2,1,:)= l;
            obj.C(2,2,:)= 2*m+l;
            obj.C(3,3,:)= m;
        end
        
    end
    
%     methods (Access = private)
%         
%         function computeYoungModulus(obj)
%             k = obj.kappa;
%             m = obj.mu;            
%             E = 4*k.*m./(k + m);
%             obj.E = E;
%         end
%         
%         function computePoissonRatio(obj)
%             k = obj.kappa;
%             m = obj.mu;            
%             nu = (k - m)./(k + m);
%             obj.nu = nu;
%         end                
%         
%     end
    
end

