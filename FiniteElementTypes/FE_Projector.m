classdef FE_Projector < Projector
    
    properties (Access = private)
        finiteElement
        FE_function
    end
    
    methods (Access = public)

        function obj = FE_Projector(cParams)
            obj.init(cParams);
            
            obj.finiteElement = FiniteElementFactory.create(cParams.feParams);       
            obj.FE_function = FE_FunctionFactory.create(obj.mesh,1,cParams.feParams);
        end

        function xFun = project(obj, x)
            LHS = obj.computeLHS();
            RHS = obj.computeRHS(x);
            xProj = LHS\RHS;
            s.mesh    = obj.mesh;
            s.fValues = xProj;
            s.polynomialOrder = obj.finiteElement.polynomialOrder;
            xFun = FE_RaviartThomasFunction(s); % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        end

    end

    methods (Access = private)
        
        function LHS = computeLHS(obj)
            s.mesh  = obj.mesh;
            s.fun = obj.FE_function;            
            s.quadratureOrder = 'CUBIC';
            s.type  = 'MassMatrixVect'; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end

        function RHS = computeRHS(obj,fun)
            quad = obj.createRHSQuadrature();
            xV = quad.posgp;
            dV = obj.mesh.computeDvolume(quad);
            
            obj.FE_function.interpolation.computeShapeDeriv(xV);
            shapes = permute(obj.FE_function.interpolation.shape,[1 3 2]);
            
            dofsGlob = obj.FE_function.computeDofConnectivity()';
            nDofs = max(max(dofsGlob));

            nGaus = quad.ngaus;
            nFlds = fun.ndimf;
            nDofel = size(shapes,1);
            
            fGa = fun.evaluate(xV);
            ndim = size(fGa,2)/nGaus;
            fGaus = zeros(ndim,nGaus,obj.mesh.nelem);
            for i = 1:ndim
                fGaus(i,:,:) = fGa(1,(i-1)*nGaus+1:i*nGaus,:);
            end
            
            locPointEdge = squeeze(obj.mesh.edges.localNodeByEdgeByElem(:,:,1));
            sides = zeros(obj.mesh.nelem,obj.mesh.edges.nEdgeByElem);
            for ielem=1:obj.mesh.nelem
                sides(ielem,:) = ones(1,obj.mesh.edges.nEdgeByElem)-2.*(locPointEdge(ielem,:)~=1:obj.mesh.edges.nEdgeByElem);
            end
            
            f     = zeros(nDofs,nFlds);
            for iField = 1:nFlds
                for igaus = 1:nGaus
                    dVg(:,1) = dV(igaus, :);
                    fG = squeeze(fGaus(:,igaus,:)); % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    for iDof = 1:nDofel
                        dofs = dofsGlob(:,iDof);
                        Ni = shapes(iDof,:,igaus);
                        int = Ni*fG.*dVg'.*sides(:,iDof)';
                        f(:,iField) = f(:,iField) + accumarray(dofs,int,[nDofs 1]);
                    end
                end
            end
            RHS = f;
        end

        function q = createRHSQuadrature(obj)
            ord = 'CUBIC';
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(ord);
        end
        
    end

end