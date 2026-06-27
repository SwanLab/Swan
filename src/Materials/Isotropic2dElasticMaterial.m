classdef Isotropic2dElasticMaterial < IsotropicElasticMaterial

    properties (Access = private)
        cohesiveMesh
    end

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
            IxI    = repmat(kronEye(N),[1 1 1 1 nGauss nElem]);
            C = 2*mu.*I + lambda.*IxI;
            
            if not(isempty(obj.cohesiveMesh))
                C(:,:,:,:,:,obj.cohesiveMesh.listCohesiveElems) = 0; 
                C = obj.defineLeverMaterial(C,m.mesh,nGauss);
            end

            % C = obj.defineLeverMaterial(C,m.mesh,nGauss);


        end

        function plot(obj,mesh)
            s.mesh = mesh;
            s.projectorType = 'P1D';
            proj = Projector.create(s);
            p1fun = proj.project(obj);
            p1fun.plot();
        end

        function setCohesiveMesh(obj,cM)
           obj.cohesiveMesh = cM;
        end

        function C = defineLeverMaterial(obj,C,mesh,nGauss)
            E = 1e10;
            baricenters = mesh.computeBaricenter';
            isLeverElem = baricenters(:,2) > 3.1;
            mu = obj.computeMuFromYoungAndPoisson(E,0.3);
            k = obj.computeKappaFromYoungAndPoisson(E,0.3,2);
            lambda = obj.computeLambdaFromShearAndBulk(mu,k,2);
            IxI    = tensorprod(eye(2),eye(2));
            CLever = 2*mu.*eye4D(2) + lambda.*IxI;        
            C(:,:,:,:,:,isLeverElem) = repmat(CLever,[1 1 1 1 nGauss sum(isLeverElem)]);
        end
    end
    

end

