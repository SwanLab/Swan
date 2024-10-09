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


            nDofs    = obj.mesh.nnodeElem*obj.mesh.nelem*x.ndimf;
            nodes     = 1:nDofs;
            dofConnec = reshape(nodes,obj.mesh.nnodeElem*x.ndimf,obj.mesh.nelem)';


            coord = obj.mesh.coord;
            coordC(:,1) = coord(connec,1);
            coordC(:,2) = coord(connec,2);

             coor = zeros(nDofs,obj.mesh.ndim);
                for idim = 1:obj.mesh.ndim
                    coorI = repmat(coordC(:,idim)',x.ndimf,1);
                    coor(:,idim) = coorI(:);
                end

            ndimf = x.ndimf;
            dofCoord = coor;

            xP1D = P1DiscontinuousFunction.create(obj.mesh,dofConnec,dofCoord,ndimf);            

            if strcmp(order, 'P1')
                f = x.fValues;
                fVals = f(connec,:);
            else
                LHS  = obj.computeLHS(xP1D);
                RHS   = obj.computeRHS(xP1D,x);
                fVals = LHS\RHS;
                fVals = reshape(fVals',xP1D.ndimf,[])';                           
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
            s.quadType = 2;
            s.type = 'ShapeFunction';            
            s.mesh = xP1D.mesh;
            int        = RHSintegrator.create(s);
            test       = xP1D;
            RHS        = int.compute(fun,test);
        end

    end

end