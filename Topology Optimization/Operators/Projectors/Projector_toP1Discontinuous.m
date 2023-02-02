classdef Projector_toP1Discontinuous < Projector

    properties (Access = private)
        field
    end

    methods (Access = public)

        function obj = Projector_toP1Discontinuous(cParams)
            obj.init(cParams);
            obj.createField();
        end

        function xProj = project(obj, x)
            LHS = obj.computeLHS();
            RHS = obj.computeRHS(x);
            f = LHS\RHS;
            fVals = obj.reshapeFValues(f, x.ndimf);
            s.mesh    = obj.mesh;
            s.fValues = fVals;
            xProj = P1DiscontinuousFunction(s);
        end

    end

    methods (Access = private)

        function createField(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATIC';
            s.galerkinType       = 'DISCONTINUOUS';
            obj.field = Field(s);
        end

        function LHS = computeLHS(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.field = obj.field;
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

        function RHS = computeRHS(obj,fun)
            quad = obj.createRHSQuadrature(fun);
            xV = quad.posgp;
            dV = obj.mesh.computeDvolume(quad);
            obj.mesh.interpolation.computeShapeDeriv(xV);
            shapes = permute(obj.mesh.interpolation.shape,[1 3 2]);
            conne = obj.createDiscontinuousConnectivity();

            nGaus = quad.ngaus;
            nFlds = fun.ndimf;
            nElem = obj.mesh.nelem;
            nNode = size(conne,2);
            nDofs = nElem*nNode;

            fGaus = fun.evaluate(xV);
            f     = zeros(nDofs,nFlds);
            for iField = 1:nFlds
                for igaus = 1:nGaus
                    dVg(:,1) = dV(igaus, :);
                    fG = squeeze(fGaus(iField,igaus,:));
                    for inode = 1:nNode
                        dofs = conne(:,inode);
                        Ni = shapes(inode,igaus);
                        int = Ni.*fG.*dVg;
                        f(dofs,iField) = f(dofs,iField) + int;
                    end
                end
            end
            RHS = f;
        end

        function conn = createDiscontinuousConnectivity(obj)
            nElem   = obj.mesh.nelem;
            nNodeEl = obj.mesh.nnodeElem;
            nDofs   = nElem*nNodeEl;
            dofs = 1:nDofs;
            conn = reshape(dofs, [nNodeEl,nElem])';
        end

        function q = createRHSQuadrature(obj, fun)
            ord = obj.determineQuadratureOrder(fun);
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(ord);
        end

        function fVals = reshapeFValues(obj, x, nFlds)
            nElem = obj.mesh.nelem;
            nNode = obj.mesh.nnodeElem;
            fVals = reshape(x',nFlds, nNode, nElem);
        end

    end

end