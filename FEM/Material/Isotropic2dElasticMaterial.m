classdef Isotropic2dElasticMaterial < IsotropicElasticMaterial

    methods (Access = public)
        
        function obj = Isotropic2dElasticMaterial(cParams)
            obj.init(cParams);          
            obj.ndim = 2;
        end
        
    end    
    
    methods (Access = protected)
        

        
        function C = computeC(obj,mu,lambda)
            m = mu;
            l = lambda;
            C = zeros(obj.nstre,obj.nstre,obj.nElem);                        
            C(1,1,:)= 2*m+l;
            C(1,2,:)= l;
            C(2,1,:)= l;
            C(2,2,:)= 2*m+l;
            C(3,3,:)= m;
            obj.C = C;
        end
        
    end

end

