classdef Isotropic3dElasticMaterial < IsotropicElasticMaterial
    
    methods (Access = public)
        
        function obj = Isotropic3dElasticMaterial(cParams)
            obj.init(cParams);
        end
        
    end
    
    methods (Access = public)

        function C = evaluate(obj,xV)
            [mu,k] = obj.computeShearAndBulk(xV);
            l = obj.computeLambdaFromShearAndBulk(mu,k,obj.ndim);
            nGaus = size(xV,2);
            nElem = size(mu,3);
            nStre = 6;
            C = zeros(nStre,nStre,nGaus,nElem);
            C(1,1,:,:) = 2*mu + l;
            C(2,2,:,:) = 2*mu + l;
            C(3,3,:,:) = 2*mu + l;
            C(1,2,:,:) = l;
            C(2,1,:,:) = l;
            C(1,3,:,:) = l;
            C(3,1,:,:) = l;
            C(3,2,:,:) = l;
            C(2,3,:,:) = l;
            C(4,4,:,:) = mu;
            C(5,5,:,:) = mu;
            C(6,6,:,:) = mu;
        end
        
    end
    
end

