classdef Projector_toP1Discontinuous < Projector

    methods (Access = public)

        function obj = Projector_toP1Discontinuous(cParams)
            obj.init(cParams);
        end

        function xProj = project(obj, x)
            if strcmp(x.order, 'P1')
                % fVals = zeros(x.nDofsElem,obj.mesh.nelem);
                % fVals = obj.reshapeFValues(fVals,x.ndimf);
                connec = obj.mesh.connec;

                f = x.fValues;
                nNode  = size(connec,2);
                nDime  = size(f,2);
                nodes = reshape(connec',1,[]);
                fe = f(nodes,:)';
                fVals = reshape(fe,nDime,nNode,[]);
            else
                LHS = obj.computeLHS();
                RHS = obj.computeRHS(x);
                f = LHS\RHS;
                fVals = obj.reshapeFValues(f, x.ndimf);
            end
            s.mesh    = obj.mesh;
            s.fValues = fVals;
            xProj = P1DiscontinuousFunction(s);
        end

    end

    methods (Access = private)

        function LHS = computeLHS(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = P1DiscontinuousFunction.create(obj.mesh, 1);
            s.trial = P1DiscontinuousFunction.create(obj.mesh, 1);
            s.quadratureOrder = 'QUADRATIC';
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

        function RHS = computeRHS(obj,fun)
            quad = obj.createRHSQuadrature(fun);
            xV = quad.posgp;
            dV = obj.mesh.computeDvolume(quad);

            trial = P1DiscontinuousFunction.create(obj.mesh, 1);
            shapes = trial.computeShapeFunctions(xV);

           % shapes = permute(obj.mesh.interpolation.shape,[1 3 2]);
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