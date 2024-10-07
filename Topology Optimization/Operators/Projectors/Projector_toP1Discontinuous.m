classdef Projector_toP1Discontinuous < Projector

    methods (Access = public)

        function obj = Projector_toP1Discontinuous(cParams)
            obj.init(cParams);
        end

        function xProj = project(obj, x)
            if isprop(x,'order')
                order = x.order;
            else
                order = [];
            end
            if strcmp(order, 'P1')
                % fVals = zeros(x.nDofsElem,obj.mesh.nelem);
                % fVals = obj.reshapeFValues(fVals,x.ndimf);
                connec = obj.mesh.connec;

                f = x.fValues;
                nNode  = size(connec,2);
                nDime  = size(f,2);
                connec = reshape(connec',1,[]);
                
                nNodes = obj.mesh.nnodeElem*obj.mesh.nelem;
                nodes  = 1:nNodes;
                newConnec = reshape(nodes,obj.mesh.nnodeElem,obj.mesh.nelem)';
           
                fe = f(connec,:);
              %  fVals = reshape(fe,nDime,nNode,[]);
                coord = obj.mesh.coord;
                coordD(:,1) = coord(connec,1);
                coordD(:,2) = coord(connec,2);
                fVals = fe;
                
            else
                LHS = obj.computeLHS();
                RHS = obj.computeRHS(x);
                f = LHS\RHS;
                fVals = obj.reshapeFValues(f, size(f,2));
            end
            s.mesh     = obj.mesh;
            s.fValues  = fVals;
            s.dofCoord = coordD;
            s.dofConnec = newConnec;
            xProj = P1DiscontinuousFunction(s);
        end

    end

    methods (Access = private)

        function LHS = computeLHS(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = P1DiscontinuousFunction.create(obj.mesh, 1);
            s.trial = P1DiscontinuousFunction.create(obj.mesh, 1);
            s.quadratureOrder = 2; % ?
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
            nElem = obj.mesh.nelem;
            nNode = size(conne,2);
            nDofs = nElem*nNode;

            fGaus = fun.evaluate(xV);
            nFlds = size(fGaus,1);
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
            q = Quadrature.create(obj.mesh,ord);
        end

        function fVals = reshapeFValues(obj, x, nFlds)
            nElem = obj.mesh.nelem;
            nNode = obj.mesh.nnodeElem;
            fVals = reshape(x',nFlds, nNode, nElem);
        end

    end

end