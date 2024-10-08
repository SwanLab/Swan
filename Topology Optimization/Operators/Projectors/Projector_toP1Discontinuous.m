classdef Projector_toP1Discontinuous < Projector

    methods (Access = public)

        function obj = Projector_toP1Discontinuous(cParams)
            obj.init(cParams);
        end

        function xP1D = project(obj, x)
            if isprop(x,'order')
                order = x.order;
            else
                order = [];
            end
            
            connec = obj.mesh.connec;
            %connec = connec(:);
            connec = reshape(connec',1,[]);

            nNodes    = obj.mesh.nnodeElem*obj.mesh.nelem;
            nodes     = 1:nNodes;
            dofConnec = reshape(nodes,obj.mesh.nnodeElem,obj.mesh.nelem)';


            coord = obj.mesh.coord;
            dofCoord(:,1) = coord(connec,1);
            dofCoord(:,2) = coord(connec,2);
            ndimf = x.ndimf;

            xP1D = P1DiscontinuousFunction.create(obj.mesh,dofConnec,dofCoord,ndimf);            

            if strcmp(order, 'P1')
                f = x.fValues;
                fVals = f(connec,:);
            else
                LHS  = obj.computeLHS(xP1D);
                RHS   = obj.computeRHS(xP1D,x);
                fVals = LHS\RHS;
                fVals = reshape(fVals,[], xP1D.ndimf);                           
            end
               xP1D.fValues  = fVals;
        end

    end

    methods (Access = private)

        function LHS = computeLHS(obj,xP1D)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = xP1D;
            s.trial = xP1D.copy();
            s.quadratureOrder = 2; % ?
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

        function RHS = computeRHS(obj,xP1D,fun)
            quad = obj.createRHSQuadrature(fun);
            xV = quad.posgp;
            dV = obj.mesh.computeDvolume(quad);

            trial  = xP1D.copy();
            shapes = trial.computeShapeFunctions(xV);

           % shapes = permute(obj.mesh.interpolation.shape,[1 3 2]);
            conne = trial.getDofConnec();%obj.createDiscontinuousConnectivity();

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