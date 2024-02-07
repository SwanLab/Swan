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
            s.quadratureOrder = 'QUADRATIC'; % no
            s.type  = 'MassMatrix';
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

        function RHS = computeRHS(obj,fun)
            quad = obj.createRHSQuadrature(fun);
            xV = quad.posgp;
            dV = obj.mesh.computeDvolume(quad);
            
            f = LagrangianFunction.create(obj.mesh, 1,obj.order);
            shapes = f.computeShapeFunctions(quad);
            conne = f.computeDofConnectivity()';

            nGaus = quad.ngaus;
            nFlds = fun.ndimf;
            nNode = size(conne,2);
            nDofs = max(max(conne));

            fGaus = fun.evaluate(xV);
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
                q = fun.quadrature;
            else
                ord = obj.determineQuadratureOrder(fun);
                quadratureOrder = 'QUADRATIC'; % no
                
                q = Quadrature.set(obj.mesh.type);
                q.computeQuadrature(quadratureOrder);
            end
        end
        
    end

end