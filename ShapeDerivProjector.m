classdef ShapeDerivProjector < Projector
    
    methods (Access = public)

        function obj = ShapeDerivProjector(cParams)
            obj.init(cParams);
        end

        function xFun = project(obj,f,df)
            LHS     = obj.computeLHS();
            RHSdfun  = obj.computeRHSdFun(df);
           % for i = 1:length(df)
           %     RHSdfun(:,i) = obj.computeRHSdFun2(df{i});
           % end
            % 
            RHSfun = obj.computeRHSFun(f);
            RHS     = RHSfun + RHSdfun;
            xProj = LHS\RHS;
            s.mesh    = obj.mesh;
            s.fValues = xProj;
            xFun = P1Function(s);
        end

    end

    methods (Access = private)
        
        function LHS = computeLHS(obj)
            LHSM    = obj.computeMassMatrix();
            LHSK    = obj.computeStiffnessMatrix();
            epsilon = obj.mesh.computeMeanCellSize();
            LHS     = LHSM + epsilon^2*LHSK;
        end

        function LHSM = computeMassMatrix(obj)
            s.test  = P1Function.create(obj.mesh,1);
            s.trial = P1Function.create(obj.mesh,1);
            s.mesh  = obj.mesh;
            s.type  = 'MassMatrix';
            s.quadratureOrder = 'QUADRATIC';
            lhs = LHSintegrator.create(s);
            LHSM = lhs.compute();
        end

        function LHSK = computeStiffnessMatrix(obj)
            s.type = 'StiffnessMatrix';
            s.mesh = obj.mesh;
            s.test  = P1Function.create(obj.mesh, 1);
            s.trial = P1Function.create(obj.mesh, 1);
            s.quadratureOrder = 'CONSTANT';
            lhs = LHSintegrator.create(s);
            LHSK = lhs.compute();
        end

        function RHS = computeRHSFun(obj,fun)
            s.type     = 'ShapeDivergence';
            s.mesh     = obj.mesh;
            s.quadratureOrder = 'QUADRATIC';            
            rhs = RHSintegrator.create(s);
            rhs = rhs.compute(fun);
            RHS = rhs.fValues;
        end

        function RHS = computeRHSdFun2(obj,dfun)
            s.type     = 'ShapeFunction';
            s.mesh     = obj.mesh;
            s.quadratureOrder = 'QUADRATIC';            
            rhs = RHSintegrator.create(s);
            RHS = rhs.compute(dfun);
        end        

        function RHS = computeRHSdFun(obj,dfun)
            quad = obj.createRHSQuadrature(dfun);
            xV = quad.posgp;
            dV = obj.mesh.computeDvolume(quad);
            obj.mesh.interpolation.computeShapeDeriv(xV);
            shapes = permute(obj.mesh.interpolation.shape,[1 3 2]); 
            conne = obj.mesh.connec;

            nGaus = quad.ngaus;
            nFlds = dfun.ndimf;
            nNode = size(conne,2);
            nDofs = obj.mesh.nnodes;

            fGaus = dfun.evaluate(xV);
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
            ord = obj.determineQuadratureOrder(fun);
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(ord);
        end
        
    end

end