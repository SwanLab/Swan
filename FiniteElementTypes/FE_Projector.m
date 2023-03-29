classdef FE_Projector < Projector
    
    properties (Access = private)
        ls
    end
    
    methods (Access = public)

        function obj = FE_Projector(cParams)
            obj.init(cParams);
            obj.ls = LagrangeElement.create("SIMPLICIAL",cParams.polynomialOrder,2);
        end

        function xFun = project(obj, x)
            LHS = obj.computeLHS();
            RHS = obj.computeRHS(x);
            xProj = LHS\RHS;
            s.mesh    = obj.mesh;
            s.fValues = xProj;
            s.polynomialOrder = obj.ls.polynomialOrder;
            xFun = FE_LagrangianFunction(s);
        end

    end

    methods (Access = private)
        
        function LHS = computeLHS(obj)
            s.mesh  = obj.mesh;
            s.fun   = FE_LagrangianFunction.create(obj.mesh, 1,obj.ls.polynomialOrder);
            s.quadratureOrder = 'CUBIC';
            s.type  = 'MassMatrix';
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

        function RHS = computeRHS(obj,fun)
            quad = obj.createRHSQuadrature(fun);
            xV = quad.posgp;
            dV = obj.mesh.computeDvolume(quad);
            
            cParams.mesh.type = obj.mesh.type;
            cParams.order = 'LINEAR';
            cParams.polynomialOrder = obj.ls.polynomialOrder;
            I = FE_Interpolation(cParams);
            I.computeShapeDeriv(xV);
            shapes = permute(I.shape,[1 3 2]);
            
            FF   = FE_LagrangianFunction.create(obj.mesh, 1,I.lagrangeElement.polynomialOrder);
            dofsGlob = FF.computeDofConnectivity()';
            nDofs = max(max(dofsGlob));

            nGaus = quad.ngaus;
            nFlds = fun.ndimf;
            nDofel = size(shapes,1);

            fGaus = fun.evaluate(xV);
            f     = zeros(nDofs,nFlds);
            for iField = 1:nFlds
                for igaus = 1:nGaus
                    dVg(:,1) = dV(igaus, :);
                    fG = squeeze(fGaus(iField,igaus,:));
                    for inode = 1:nDofel
                        dofs = dofsGlob(:,inode);
                        Ni = shapes(inode,1,igaus);
                        int = Ni*fG.*dVg;
                        f(:,iField) = f(:,iField) + accumarray(dofs,int,[nDofs 1]);
                    end
                end
            end
            RHS = f;
        end

        function q = createRHSQuadrature(obj, fun)
%             ord = obj.determineQuadratureOrder(fun);ç
            ord = 'CUBIC';
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(ord);
        end
        
    end

end