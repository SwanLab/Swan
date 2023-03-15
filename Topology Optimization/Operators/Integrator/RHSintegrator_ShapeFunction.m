classdef RHSintegrator_ShapeFunction < handle

    properties (Access = private)
        npnod
        mesh
        globalConnec
        fType
        fNodal
        fGauss
        xGauss
        quadOrder
        quadrature
    end

    methods (Access = public)

        % Via Integrator_Simple + Integrator
        function obj = RHSintegrator_ShapeFunction(cParams)
            obj.init(cParams);
        end

        function rhs = compute(obj, fNodal)
            obj.computeFgaussFromFnodal(fNodal);
            rhsElem = obj.computeElementalRHS();
            rhs = obj.assembleIntegrand(rhsElem);
        end

        function rhs = computeFromFgauss(obj, fGaus, xGaus)
            obj.fGauss = fGaus;
            obj.xGauss = xGaus;
            rhsElem = obj.computeElementalRHS();
            rhs = obj.assembleIntegrand(rhsElem);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh         = cParams.mesh;
            obj.npnod        = cParams.npnod;
            obj.quadOrder    = 'LINEAR';
            obj.globalConnec = cParams.globalConnec;
            obj.quadrature   = obj.computeQuadrature();
        end

        function q = computeQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(obj.quadOrder);
        end

        function computeFgaussFromFnodal(obj, fNodal)
            obj.computeGaussPoints();
            obj.computeFgauss(fNodal);
        end

        function computeGaussPoints(obj)
            q = obj.quadrature;
            xG = repmat(q.posgp,[1,1,obj.mesh.nelem]);
            obj.xGauss = xG;
        end

        function computeFgauss(obj, fNodal)
%             mmm.coord = obj.mesh.coord;
%             mmm.connec = obj.globalConnec;
%             locMesh = Mesh(mmm);
            msh.connec = obj.globalConnec;
            msh.type   = obj.mesh.type;
            s.mesh    = msh;
            s.fValues = fNodal;
%             s.mesh    = locMesh;
%             s.connec = obj.globalConnec;
%             s.type   = obj.mesh.type;
            f = P1Function(s);
            fG = f.evaluate(obj.xGauss);
            fG = permute(fG,[2 3 1]);
            obj.fGauss = fG;
        end

        function rhsC = computeElementalRHS(obj) % integrate@RHSintegrator
            fG     = obj.fGauss;
            dV     = obj.computeDvolume();
%             fdV    = (fG.*dV);
            shapes = obj.computeShapeFunctions();
            nnode  = size(shapes,1);
            nelem  = size(shapes,2);
            int = zeros(nnode,nelem);
            for igaus = 1:obj.quadrature.ngaus
                fdv = fG(igaus,:).*dV(igaus,:);
                shape = shapes(:, :, igaus);
                int = int + bsxfun(@times,shape,fdv);
            end
            rhsC = transpose(int);
        end

        function f = assembleIntegrand(obj,rhsElem)
            integrand = rhsElem;
            ndofs = obj.npnod;
            connec = obj.globalConnec;
            nnode  = size(connec,2);
            f = zeros(ndofs,1);
            for inode = 1:nnode
                int = integrand(:,inode);
                con = connec(:,inode);
                f = f + accumarray(con,int,[ndofs,1],@sum,0);
            end
        end

        function dV = computeDvolume(obj)
            q = obj.quadrature;
            dV = obj.mesh.computeDvolume(q);
        end

        function shapes = computeShapeFunctions(obj)
            int = Interpolation.create(obj.mesh,'LINEAR');
            int.computeShapeDeriv(obj.xGauss);
            shapes = permute(int.shape,[1 3 2]);
        end

    end

end
