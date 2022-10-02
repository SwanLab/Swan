classdef Projector_toP1Discontinuous < Projector

    properties (Access = public)

    end

    properties (Access = private)
        mesh
        connec
        quadOrder
    end

    properties (Access = private)
        meshD
        field
    end

    methods (Access = public)

        function obj = Projector_toP1Discontinuous(cParams)
            obj.init(cParams);
            obj.createDiscontinuousMesh();
            obj.createField();
        end

        function xProj = project(obj, x)
            LHS = obj.computeLHS();
            RHS = obj.computeRHS(x);
            f = LHS\RHS;
            fVals = obj.reshapeFValues(f, x.ndimf);
            s.type    = obj.mesh.type;
            s.connec  = obj.mesh.connec;
            s.fValues = fVals;
            xProj = P1DiscontinuousFunction(s);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh   = cParams.mesh;
            obj.connec = cParams.connec;
            obj.quadOrder = 'QUADRATIC'; % Should change according to fun
        end

        function createDiscontinuousMesh(obj)
            obj.meshD = obj.mesh.createDiscontinuousMesh();
        end

        function createField(obj)
            s.mesh               = obj.meshD;
            s.ndimf              = 1;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATIC';
            s.
            obj.field = Field(s);
        end

        function LHS = computeLHS(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.meshD;
            s.field = obj.field;
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();

            %             sK.type = 'StiffnessMatrix';
            %             sK.mesh  = obj.meshD;
            %             sK.field = obj.field;
            %             lhsK = LHSintegrator.create(sK);
            %             LHSK = lhsK.compute();
            %
            %             epsilon = obj.mesh.computeMeanCellSize();
            %             LHS = LHS + LHSK*epsilon^2;
        end

        function RHS = computeRHS(obj,fun)
            quad = obj.createRHSQuadrature(fun);
            xV = quad.posgp;
            dV = obj.mesh.computeDvolume(quad);
            obj.mesh.interpolation.computeShapeDeriv(xV);
            shapes = permute(obj.mesh.interpolation.shape,[1 3 2]);
            conne = obj.meshD.connec;

            nGaus = quad.ngaus;
            nFlds = fun.ndimf;
            nElem = obj.mesh.nelem;
            nNode = size(conne,2);
            nDofs = obj.meshD.nnodes;

            fLoc = zeros(nNode,nElem,nFlds);
            fGaus = fun.evaluate(xV);
            %             for iElem = 1:nElem
            for igaus = 1:nGaus
                dVg(:,1) = dV(igaus, :);
                for iField = 1:nFlds
                    fG = squeeze(fGaus(iField,igaus,:));
                    for inode = 1:nNode
                        Ni = shapes(inode,igaus);
                        f(:,1) = squeeze(fLoc(inode,:,iField));
                        int = Ni.*fG.*dVg;
                        fLoc(inode,:,iField) = f + int;
                    end
                end
            end
            %             end
            f = zeros(nDofs,nFlds);
            for iField = 1:nFlds
                for iElem = 1:nElem
                    for inode = 1:nNode
                        dofs = conne(iElem,inode);
                        f(dofs,iField) = fLoc(inode,iElem,iField);
                    end
                end
            end
            RHS = f;
        end

        function q = createRHSQuadrature(obj, fun)
            ord = obj.determineQuadratureOrder(fun);
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(ord);
        end

        function fVals = reshapeFValues(obj, x, nFlds)
            nElem = obj.meshD.nelem;
            nNode = obj.meshD.nnodeElem;
            %             fRshp = reshape(x, [nNode, nFlds, nElem]);
            %             fVals = permute(fRshp, [2 1 3]);
            fVals = reshape(x',nFlds, nNode, nElem);
        end

    end

end

