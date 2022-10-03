classdef H1Projector_toP1Discontinuous < Projector

    properties (Access = private)
        fieldMass
        fieldStiffness
    end

    methods (Access = public)

        function obj = H1Projector_toP1Discontinuous(cParams)
            obj.init(cParams);
            obj.createFieldMass();
            obj.createFieldStiffness();
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

        function createFieldMass(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATIC';
            s.galerkinType       = 'DISCONTINUOUS';
            obj.fieldMass = Field(s);
        end

        function createFieldStiffness(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'CONSTANT';
            s.galerkinType       = 'DISCONTINUOUS';
            obj.fieldStiffness = Field(s);
        end

        function LHS = computeLHS(obj)
            LHSM    = obj.computeMassMatrix();
            LHSK    = obj.computeStiffnessMatrix();
            epsilon = obj.mesh.computeMeanCellSize();
            LHS     = LHSM + epsilon^2*LHSK;
        end

        function LHSM = computeMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.field = obj.fieldMass;
            lhs = LHSintegrator.create(s);
            LHSM = lhs.compute();
        end

        function LHSK = computeStiffnessMatrix(obj)
            s.type  = 'StiffnessMatrix';
            s.mesh  = obj.mesh;
            s.field = obj.fieldStiffness;
            lhs = LHSintegrator.create(s);
            LHSK = lhs.compute();
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

            fLoc  = zeros(nNode,nElem,nFlds);
            fGaus = fun.evaluate(xV);
            f     = zeros(nDofs,nFlds);
            for iField = 1:nFlds
                for igaus = 1:nGaus
                    dVg(:,1) = dV(igaus, :);
                    fG = squeeze(fGaus(iField,igaus,:));
                    for inode = 1:nNode
                        Ni = shapes(inode,igaus);
                        fL(:,1) = squeeze(fLoc(inode,:,iField));
                        int = Ni.*fG.*dVg;
                        fLoc(inode,:,iField) = fL + int;
                    end
                end
                for iElem = 1:nElem
                    dofs = conne(iElem,:);
                    f(dofs,iField) = fLoc(:,iElem,iField);
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