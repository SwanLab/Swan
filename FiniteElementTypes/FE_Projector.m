classdef FE_Projector < Projector
    
    properties (Access = private)
        finiteElement
        FE_function
        elementType
    end
    
    methods (Access = public)

        function obj = FE_Projector(cParams)
            obj.init(cParams);
            
            obj.finiteElement = FiniteElementFactory.create(cParams.feParams);
            cParams.feParams.finiteElement = obj.finiteElement;
            obj.FE_function = FE_FunctionFactory.create(obj.mesh,1,cParams.feParams);
            obj.elementType = cParams.feParams.type;
        end

        function xFun = project(obj, x)
            LHS = obj.computeLHS();
            RHS = obj.computeRHS(x);
            xProj = LHS\RHS;
            s.mesh    = obj.mesh;
            s.fValues = xProj;
            s.polynomialOrder = obj.finiteElement.polynomialOrder;
            s.type = obj.elementType;
            s.finiteElement = obj.finiteElement;
            xFun = FE_FunctionFactory.createWithFValues(s);
        end

    end

    methods (Access = private)
        
        function LHS = computeLHS(obj)
            s.mesh  = obj.mesh;
            s.fun = obj.FE_function;            
            s.quadratureOrder = 'ORDER4';
            s.type  = 'MassMatrix'; % For Lagrange
%             s.type  = 'MassMatrixVect'; % For Raviart-Thomas or Nedelec
%             s.type  = 'CurlNCurlN';
%             s.type  = 'ElasticStiffnessMatrix';
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
        end
        
        function RHS = computeRHS(obj,fun)
            switch obj.elementType
                case "Lagrange Simplicial"
                    RHS = obj.computeRHS_LS(fun);
                case "Raviart-Thomas"
                    RHS = obj.computeRHS_RT(fun);
                case "Nedelec"
                    RHS = obj.computeRHS_N(fun);
            end
        end
        
        function RHS = computeRHS_LS(obj,fun)
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
            fGaus = fun.evaluate(xV);
            
            f     = zeros(nDofs,nFlds);
            for iField = 1:nFlds
                for igaus = 1:nGaus
                    dVg(:,1) = dV(igaus, :);
                    fG = squeeze(fGaus(iField,igaus,:));
                    for iDof = 1:nDofel
                        dofs = dofsGlob(:,iDof);
                        Ni = shapes(iDof,1,igaus);
                        int = Ni*fG.*dVg;
                        f(:,iField) = f(:,iField) + accumarray(dofs,int,[nDofs 1]);
                    end
                end
            end
            RHS = f;
        end
        
        function RHS = computeRHS_RT(obj,fun)
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
            
            J = obj.mesh.geometry.jacobian;
            Jdet = obj.mesh.geometry.detJ;
            
            f     = zeros(nDofs,nFlds);
            for iField = 1:nFlds
                for igaus = 1:nGaus
                    dVg(:,1) = dV(igaus, :);
                    fG = squeeze(fGaus(:,igaus,:));
                    for iDof = 1:nDofel
                        Jd = Jdet(:,igaus);
                        dofs = dofsGlob(:,iDof);
                        Ni = pagemtimes(squeeze(shapes(iDof,:,igaus)),J);
                        if (size(squeeze(Ni),1)~=1) 
                            Ni = squeeze(Ni)';
                        end
                        int = diag(Ni*fG)'.*(dVg./Jd)'.*sides(:,iDof)';
                        f(:,iField) = f(:,iField) + accumarray(dofs,int,[nDofs 1]);
                    end
                end
            end
            RHS = f;
        end

        function RHS = computeRHS_N(obj,fun)
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
            
            R = [0,1;-1,0];
            J = obj.mesh.geometry.jacobian;
            J = pagemtimes(pagemtimes(R,J),R');
            Jdet = obj.mesh.geometry.detJ;
            
            f     = zeros(nDofs,nFlds);
            for iField = 1:nFlds
                for igaus = 1:nGaus
                    dVg(:,1) = dV(igaus, :);
                    fG = squeeze(fGaus(:,igaus,:)); % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    for iDof = 1:nDofel
                        Jd = Jdet(:,igaus);
                        dofs = dofsGlob(:,iDof);
                        Ni = pagemtimes(squeeze(shapes(iDof,:,igaus)),J);
                        if (size(squeeze(Ni),1)~=1) 
                            Ni = squeeze(Ni)';
                        end
                        int = diag(Ni*fG)'.*(dVg./Jd)'.*sides(:,iDof)';
                        f(:,iField) = f(:,iField) + accumarray(dofs,int,[nDofs 1]);
                    end
                end
            end
            RHS = f;
        end

        function q = createRHSQuadrature(obj)
            ord = 'ORDER10';
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(ord);
        end
        
    end

end