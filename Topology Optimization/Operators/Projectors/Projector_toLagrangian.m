classdef Projector_toLagrangian < Projector
    
    properties (Access = private)
        order
    end
    
    methods (Access = public)

        function obj = Projector_toLagrangian(cParams)
            obj.init(cParams);
            obj.order = cParams.projectorType;
        end

        function xFun = project(obj, x)
            LHS = obj.computeLHS();
            RHS = obj.computeRHS(x);
            xProj = LHS\RHS;
            s.mesh    = obj.mesh;
            s.fValues = xProj;
            s.order = obj.order;
            xFun = LagrangianFunction(s);
        end

    end

    methods (Access = private)
        
        function LHS = computeLHS(obj)
            s.mesh  = obj.mesh;
            s.test  = LagrangianFunction.create(obj.mesh, 1, obj.order);
            s.trial = LagrangianFunction.create(obj.mesh, 1, obj.order);
            s.quadratureOrder = 'ORDER10'; % no
            s.type  = 'MassMatrix';
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

        function RHS = computeRHS(obj,fun)
            quad = obj.createRHSQuadrature(fun);
            xV = quad.posgp;
            fGaus = fun.evaluate(xV);
            dV = obj.mesh.computeDvolume(quad);
            
            f = LagrangianFunction.create(obj.mesh, 1,obj.order);
            shapes = f.computeShapeFunctions(xV);
            conne = f.computeDofConnectivity()';

            nGaus = size(xV,2);
            nFlds = size(fGaus,1);
            nNode = size(conne,2);
            nDofs = max(max(conne));

            f     = zeros(nDofs,nFlds);
            for iField = 1:nFlds
                for igaus = 1:nGaus
                    dVg(:,1) = dV(igaus, :);
                    fG = squeeze(fGaus(iField,igaus,:));
                    for inode = 1:nNode
                        dofs = conne(:,inode);
                        Ni = shapes(inode,igaus);
                        int = Ni*fG.*dVg;
                        f(:,iField) = f(:,iField) + accumarray(dofs,int,[nDofs 1]);
                    end
                end
            end
            RHS = f;
        end

        function q = createRHSQuadrature(obj, fun)
            if isa(fun, 'FGaussDiscontinuousFunction')
                ord = fun.getQuadratureOrder();
            else
                ord = obj.determineQuadratureOrder(fun);
                ord = 'ORDER10'; % no
            end
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(ord);
        end

    end

end