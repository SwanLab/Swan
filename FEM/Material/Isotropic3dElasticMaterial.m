classdef Isotropic3dElasticMaterial < IsotropicElasticMaterial
    
    methods (Access = public)
        
        function obj = Isotropic3dElasticMaterial(cParams)
            obj.init(cParams);
        end
        
    end
    
    methods (Access = public)

        function C = compute(obj)
            m = obj.shear;
            l = obj.lambda;     
            nGaus = size(m,1);                        
            nElem = size(m,2);
            C = zeros(obj.nstre,obj.nstre,nGaus,nElem);
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
        end
        
    end
    
end

