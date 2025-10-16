classdef Isotropic2dElasticMaterial < IsotropicElasticMaterial

    methods (Access = public)
        
        function obj = Isotropic2dElasticMaterial(cParams)
            obj.init(cParams);
        end

        function C = evaluate(obj,xV)
            [l,m] = computeLameParameters(obj);
            lambda = l.evaluate(xV);
            mu = m.evaluate(xV);

            N = obj.ndim;
            nGauss = size(mu,2);
            nElem  = m.mesh.nelem;
            lambda = reshape(lambda,[1 1 1 1 nGauss nElem]);
            mu     = reshape(mu,[1 1 1 1 nGauss nElem]);
            I      = repmat(eye4D(N),[1 1 1 1 nGauss nElem]);
            IxI    = permute(repmat(kronEye(N),[1 1 1 1 nGauss nElem]),[4 3 2 1 5 6]);
            C = 2*mu.*I + lambda.*IxI;
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

