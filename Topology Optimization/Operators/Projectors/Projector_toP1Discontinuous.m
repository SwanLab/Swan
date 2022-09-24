classdef Projector_toP1Discontinuous < handle

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
        quadrature
    end

    methods (Access = public)

        function obj = Projector_toP1Discontinuous(cParams)
            obj.init(cParams);
            obj.createDiscontinuousMesh();
            obj.createQuadrature();
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
            obj.quadOrder = 'LINEAR';
        end

        function createDiscontinuousMesh(obj)
            obj.meshD = obj.mesh.createDiscontinuousMesh();
        end

        function q = createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(obj.quadOrder);
            obj.quadrature = q;
        end

        function createField(obj)
            s.mesh               = obj.meshD;
            s.ndimf              = 1;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATIC';
            obj.field = Field(s);
        end

        function LHS = computeLHS(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.meshD;
            s.field = obj.field;
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

        function RHS = computeRHS(obj,fun)
            xV = obj.quadrature.posgp;
            dV = obj.meshD.computeDvolume(obj.quadrature);
            obj.meshD.interpolation.computeShapeDeriv(xV);
            shapes = permute(obj.meshD.interpolation.shape,[1 3 2]);
            conne = obj.meshD.connec;

            nGaus = obj.quadrature.ngaus;
            nFlds = fun.ndimf;
            nElem = obj.meshD.nelem;
            nNods = size(shapes,1);
            nNode = size(conne,2);
            nDofs = obj.meshD.nnodes;

            rhs = zeros(nNods,nElem, nFlds);
            f = zeros(nDofs,nFlds);
            fGaus = fun.evaluate(xV);
            for igaus = 1:nGaus
                dVg(:,1) = dV(igaus, :);
                for iField = 1:nFlds
                    fG = squeeze(fGaus(iField,:,:));
                    fdVg = fG.*dVg;
                    Ni = shapes(:,:, igaus);
                    rhs(:,:,iField) = rhs(:,:,iField) + bsxfun(@times,Ni,fdVg');
                    for inode = 1:nNode
                        int = rhs(inode,:,iField);
                        con = conne(:,inode);
                        f(:,iField) = f(:,iField) + accumarray(con,int,[nDofs,1],@sum,0);
                    end
                end
            end
            RHS = f;
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

