classdef Isotropic2dElasticMaterial < IsotropicElasticMaterial

    methods (Access = public)
        
        function obj = Isotropic2dElasticMaterial(cParams)
            obj.init(cParams);
        end

        function C = evaluate(obj,xV)
            [mu,k] = obj.computeShearAndBulk(xV);
            l = obj.computeLambdaFromShearAndBulk(mu,k,obj.ndim);
            nGaus = size(xV,2);
            nElem = size(mu,3);
            nStre = 3;
            C = zeros(nStre,nStre,nGaus,nElem);
            C(1,1,:,:)= 2*mu+l;
            C(1,2,:,:)= l;
            C(2,1,:,:)= l;
            C(2,2,:,:)= 2*mu+l;
            C(3,3,:,:)= mu;
        end

        function plot(obj,mesh)
            s.mesh = mesh;
            s.projectorType = 'P1D';
            proj = Projector.create(s);
            p1fun = proj.project(obj);
            p1fun.plot();
        end
        
    end
    

end

