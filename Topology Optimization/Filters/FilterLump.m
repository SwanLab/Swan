classdef FilterLump < handle

    properties (Access = private)
        trial
        mesh
    end

    methods (Access = public)

        function obj = FilterLump(cParams)
            obj.init(cParams);
        end

        function xFun = compute(obj, x, quadType)
            LHS = obj.computeLHS();
            LHS = diag(sum(LHS));
            RHS = obj.computeRHS(x, quadType);
            xProj = LHS\RHS;
            s.mesh    = obj.mesh;
            s.fValues = xProj;
            xFun = P1Function(s);
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            cParams.feFunType = 'P1Function'; % class(trial) will come from outside
            cParams.ndimf     = 1;
            obj.trial         = FeFunction.createEmpty(cParams);
            obj.mesh          = cParams.mesh;
        end

        function LHS = computeLHS(obj)
            s.mesh  = obj.mesh;
            s.test  = P1Function.create(obj.mesh, 1);
            s.trial = P1Function.create(obj.mesh, 1);
            s.quadratureOrder = 'QUADRATIC';
            s.type  = 'MassMatrix';
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

        function RHS = computeRHS(obj,fun, quadType)
            quad = obj.createRHSQuadrature(quadType);
            xV = quad.posgp;
            dV = obj.mesh.computeDvolume(quad);
            obj.mesh.interpolation.computeShapeDeriv(xV);
            
            shapes = obj.trial.computeShapeFunctions(quad);

            conne = obj.mesh.connec;

            nGaus = quad.ngaus;
            nFlds = fun.ndimf;
            nNode = size(conne,2);
            nDofs = obj.mesh.nnodes;

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

        function q = createRHSQuadrature(obj, quadType)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(quadType);
        end
        
    end

end