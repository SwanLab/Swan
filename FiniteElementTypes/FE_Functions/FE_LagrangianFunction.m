classdef FE_LagrangianFunction < FeFunction

    properties (Access = public)
        dofs
        ndofs
        coords
        polynomialOrder
        interpolation
    end

    properties (Access = private)
        meshCoarse
        meshFine
    end

    methods (Access = public)

        function obj = FE_LagrangianFunction(cParams)
            obj.init(cParams);
            obj.createInterpolation();
            obj.computeDofs();
            obj.computeCoords();
        end

        function fxV = evaluate(obj, xV)
            obj.interpolation.computeShapeDeriv(xV);
            shapes = obj.interpolation.shape;
            nNode  = size(shapes,1);
            nGaus  = size(shapes,2);
            nF     = size(obj.fValues,2);
            nElem  = size(obj.mesh.connec,1);
            fxV = zeros(nF,nGaus,nElem);
            for iGaus = 1:nGaus
                for iNode = 1:nNode
                    node = obj.mesh.connec(:,iNode);
                    Ni = shapes(iNode,iGaus);
                    fi = obj.fValues(node,:);
                    f(:,1,:) = Ni*fi';
                    fxV(:,iGaus,:) = fxV(:,iGaus,:) + f;
                end
            end
        end

        function N = computeShapeFunctions(obj, quad)
%             obj.mesh.computeInverseJacobian(quad,obj.interpolation);
            xV = quad.posgp;
            obj.interpolation.computeShapeDeriv(xV);
            N = obj.interpolation.shape;
        end
        
        function dNdx  = computeCartesianDerivatives(obj,quad)
            switch obj.mesh.type
                case 'LINE'
                    invJ  = obj.mesh.computeInverseJacobian(quad,obj.interpolation);
                    nElem = obj.mesh.nelem;
                    nNode = obj.interpolation.nnode;
                    nDime = obj.mesh.ndim;
                    nGaus = quad.ngaus;
                    deriv  = obj.mesh.interpolation.deriv(1,:,:,:);
                    dShapes = deriv;
                    dN = zeros(nDime,nNode,nElem,nGaus);
                    for iGaus = 1:nGaus
                        invJg = invJ(:,1);
                        for iDim = 1:nDime
                            for iNode = 1:nNode
                                dShapeI(:,1) = squeeze(dShapes(1,iNode,iGaus,:));
                                dN(iDim,iNode,:,iGaus) = invJg.*dShapeI;
                            end
                        end
                    end
                    dNdx = dN;
                    
                otherwise
                    nElem = size(obj.mesh.connec,1);
                    nNode = obj.interpolation.nnode;
                    nDime = obj.interpolation.ndime;
                    nGaus = quad.ngaus;
                    invJ  = obj.mesh.computeInverseJacobian(quad,obj.interpolation);
                    dShapeDx  = zeros(nDime,nNode,nElem,nGaus);
                    for igaus = 1:nGaus
                        dShapes = obj.interpolation.deriv(:,:,igaus);
                        for jDime = 1:nDime
                            invJ_JI   = invJ(:,jDime,:,igaus);
                            dShape_KJ = dShapes(jDime,:);
                            dSDx_KI   = bsxfun(@times, invJ_JI,dShape_KJ);
                            dShapeDx(:,:,:,igaus) = dShapeDx(:,:,:,igaus) + dSDx_KI;
                        end
                    end
                    dNdx = dShapeDx;
            end
        end

        function gradFun = computeGradient(obj, quad)
            dNdx = obj.computeCartesianDerivatives(quad);
            nDimf = obj.ndimf;
            nDims = size(dNdx, 1); % derivX, derivY (mesh-related?)
            nNode = size(dNdx, 2);
            nElem = size(dNdx, 3);
            nGaus = size(dNdx, 4);
            
            grad = zeros(nDims,nDimf, nElem, nGaus);
            for iGaus = 1:nGaus
                dNdx_g = dNdx(:,:,:,iGaus);
                for iDims = 1:nDims
                    for iNode = 1:nNode
                        dNdx_i = squeeze(dNdx_g(iDims, iNode,:));
                        nodes  = obj.mesh.connec(:,iNode);
                        f      = obj.fValues(nodes,:);
                        dfdx(1,:,:) = (dNdx_i.*f)';
                        grad(iDims,:,:,iGaus) = grad(iDims,:,:,iGaus) + dfdx;
                    end
                end
            end
            fVR = reshape(grad, [nDims*nDimf,nElem, nGaus]);
            s.fValues    = permute(fVR, [1 3 2]);
            s.mesh       = obj.mesh;
            s.quadrature = quad;
            gradFun = FGaussDiscontinuousFunction(s);
        end

        function symGradFun = computeSymmetricGradient(obj,quad)
            grad = obj.computeGradient(quad);
            nDimf = obj.ndimf;
            nDims = size(grad.fValues, 1)/nDimf;
            nGaus = size(grad.fValues, 2);
            nElem = size(grad.fValues, 3);

            gradReshp = reshape(grad.fValues, [nDims,nDimf,nGaus,nElem]);
            gradT = permute(gradReshp, [2 1 3 4]);
            symGrad = 0.5*(gradReshp + gradT);
            
            s.fValues    = reshape(symGrad, [nDims*nDimf,nGaus,nElem]);
            s.quadrature = quad;
            s.mesh       = obj.mesh;
            symGradFun = FGaussDiscontinuousFunction(s);
        end

        function divF = computeDivergence(obj,q)
            dNdx = obj.computeCartesianDerivatives(q);
            fV = obj.fValues;
            nodes = obj.mesh.connec;
            nNode = obj.mesh.nnodeElem;
            nDim  = obj.mesh.ndim;
            divV = zeros(q.ngaus,obj.mesh.nelem);
            for igaus = 1:q.ngaus
                for kNode = 1:nNode
                    nodeK = nodes(:,kNode);
                    for rDim = 1:nDim
                        dNkr = squeeze(dNdx(rDim,kNode,:,igaus));
                        fkr = fV(nodeK,rDim);
                        int(1,:) = dNkr.*fkr;
                        divV(igaus,:) = divV(igaus,:) + int;
                    end
                end
            end
            s.quadrature = q;
            s.mesh       = obj.mesh;
            s.fValues(1,:,:) = divV;
            divF = FGaussDiscontinuousFunction(s);
        end

        function fdivF = computeFieldTimesDivergence(obj,q)
            fG  = obj.evaluate(q.posgp);
            dfG = obj.computeDivergence(q);
            fdivFG = bsxfun(@times,dfG.fValues,fG);
            s.quadrature = q;
            s.mesh       = obj.mesh;
            s.fValues    = fdivFG;
            fdivF = FGaussDiscontinuousFunction(s);
        end

        function fFine = refine(obj,m,mFine)
            fNodes  = obj.fValues;
            fEdges  = obj.computeFunctionInEdges(m, fNodes);
            fAll    = [fNodes;fEdges];
            s.mesh    = mFine;
            s.fValues = fAll;
            fFine = FE_LagrangianFunction(s);
        end

        function dofConnec = computeDofConnectivity(obj)
            dofConnec = obj.dofs';
        end

        function plot(obj) % 2D domains only

            switch obj.mesh.type
                case {'TRIANGLE','QUAD'}
                    figure()
                    hold on
                    for idim = 1:obj.ndimf
                        subplot(1,obj.ndimf,idim);
                        nelem = size(obj.dofs,1);
                        for ielem = 1:nelem
                            sM.coord = obj.mesh.coord(obj.mesh.connec(ielem,:),:);
                            sM.connec = [1 2 3];
                            m = Mesh(sM);
                            m = m.remesh(2); m = m.remesh(2); m = m.remesh(2);
                            x = m.coord(:,1);
                            y = m.coord(:,2);
                            z = obj.interpolation.lagrangeElement.evaluate(obj.fValues(obj.dofs(ielem,:)),x,y);
                            a = trisurf(m.connec,x,y,z); 
                            a.EdgeColor = [0 0 0];
                            a.EdgeAlpha = 0.3;
                        end
%                         shading interp
                        grid on
                        view(0,90)
                        obj.mesh.plot();
                        title(['dim = ', num2str(idim)]);
                    end
                  
                case 'LINE'
                    x = obj.mesh.coord(:,1);
                    y = obj.fValues;
                    figure()
                    plot(x,y)
            end
        end

        function plotArrowVector(obj)
            figure()
            a = obj.fValues;
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            ax = a(:,1);
            ay = a(:,2);
            q = quiver(x,y,ax,ay);
            q.ShowArrowHead = 'off';
        end

        function print(obj, s)
%             s.mesh
            s.mesh = obj.mesh;
            s.fun = {obj};
            p = FunctionPrinter(s);
            p.print();
        end

        function [res, pformat] = getDataToPrint(obj)
            nNods = size(obj.fValues, 1);
            s.nDimf   = obj.ndimf;
            s.nData   = nNods;
            s.nGroup  = nNods;
            s.fValues = obj.getFormattedFValues();
            fps = FunctionPrintingSettings(s);
            [res, pformat] = fps.getDataToPrint();
        end

    end

    methods (Access = public, Static)

        function p = create(mesh, ndimf, polOrder)
            d = FE_LagrangianFunction.numberDofs(mesh,polOrder);
            s.fValues = zeros(d, ndimf);
            s.mesh    = mesh;
            s.polynomialOrder = polOrder;
            p = FE_LagrangianFunction(s);
        end

        function fS = times(f1,f2)
            fS = f1.fValues.*f2.fValues;
            s.fValues = fS;
            s.mesh    = f1.mesh;
            fS = FE_LagrangianFunction(s);
        end
        
        function ndofs = numberDofs(mesh,order)
            if order == 1
                ndofs = mesh.nnodes;
            else
                mesh.computeEdges();
                if order >= 3
                    nodesel = nchoosek(order-1,order-3);
                else
                    nodesel = 0;
                end
                ndofs = mesh.nnodes + max(max(mesh.edges.edgesInElem))*(order-1) + mesh.nelem*nodesel;
            end
        end
        
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
            obj.fValues = cParams.fValues;
            obj.ndimf   = size(cParams.fValues,2);
            obj.order   = 'LINEAR';
            obj.meshCoarse = cParams.mesh;
            obj.polynomialOrder = cParams.polynomialOrder;
        end
        
        function [x,y,z,connec] = createPlotMesh(obj)
            nelem = size(obj.dofs,1);
            ndofsel = size(obj.dofs,2);
            n = 15;
            nn = 16;
            x = zeros(nelem*n,1);
            y = zeros(nelem*n,1);
            z = zeros(nelem*n,1);
            connec = zeros(nelem*nn,3);
            for ielem = 1:nelem
                sM.coord = obj.mesh.coord(obj.mesh.connec(ielem,:),:);
                sM.connec = [1 2 3];
                m = Mesh(sM);
                m = m.remesh(2); m = m.remesh(2);
                x((ielem-1)*n+1:ielem*n) = m.coord(:,1);
                y((ielem-1)*n+1:ielem*n) = m.coord(:,2);
                z((ielem-1)*n+1:ielem*n) = obj.interpolation.lagrangeElement.evaluate(obj.fValues(obj.dofs(ielem,:)),m.coord(:,1),m.coord(:,2));
                connec((ielem-1)*nn+1:ielem*nn,:) = m.connec+(ielem-1)*n;
            end
            
            
%             order = obj.polynomialOrder;
%             nlocel = order^2;
%             nelem = obj.mesh.nelem;
%             connec = zeros(nelem*nlocel,3);
%             for ielem = 1:size(obj.dofs,1)
%                 con = zeros(nlocel,3);
%                 switch obj.polynomialOrder
%                     case 1
%                         con = [1 2 3];
%                     case 2
%                         con  = [1 4 6;2 4 5;4 5 6;3 5 6];
%                     case 3
%                         con = [1 4 9;4 9 10;4 5 10;5 7 10;5 2 7;8 9 10;6 8 10;6 7 10;6 8 3];
%                     case 4
%                         con = [];
%                 end
%                 for i = 1:nlocel
%                     connec((ielem-1)*nlocel+i,:) = obj.dofs(ielem,con(i,:));
%                 end
%             end
        end

        function createInterpolation(obj)
            cParams.mesh.type = obj.mesh.type;
            cParams.order = obj.order;
            cParams.polynomialOrder = obj.polynomialOrder;
            obj.interpolation = FE_Interpolation(cParams);
        end

        % Printing
        function fM = getFormattedFValues(obj)
            fM = obj.fValues;
        end

        function f = computeFunctionInEdges(~,m,fNodes)
            s.edgeMesh = m.computeEdgeMesh();
            s.fNodes   = fNodes;
            eF         = EdgeFunctionInterpolator(s);
            f = eF.compute();
        end
        
        function computeDofs(obj)
            dofsVertices = obj.computeDofsVertices();
            dofsEdges = obj.computeDofsEdges();
            dofsElements = obj.computeDofsElements(dofsEdges);
            
            obj.dofs = [dofsVertices,dofsEdges,dofsElements];
            obj.ndofs = max(max(obj.dofs));
        end
        
        function dofsVertices = computeDofsVertices(obj)
            mesh = obj.mesh;
            dofsVertices = mesh.connec;
        end
        
        function dofsEdges = computeDofsEdges(obj)
            mesh = obj.mesh;
            mesh.computeEdges();
            if obj.polynomialOrder == 1
                dofsEdges = [];
            else
                edges = mesh.edges.edgesInElem;
                ndofEd = (obj.polynomialOrder-1);
                ndofsEdgeElem = ndofEd*obj.mesh.nnodeElem;
                dofsEdges = zeros(obj.mesh.nelem,ndofsEdgeElem);
                for i = 1:obj.mesh.nelem
                    for j = 1:obj.mesh.nnodeElem
                        dofsEdges(i,(j-1)*ndofEd+1:j*ndofEd) = edges(i,j)*ndofEd-(ndofEd-1):edges(i,j)*ndofEd;
                    end
                end
                dofsEdges = dofsEdges + mesh.nnodes;
            end
        end
        
        function dofsElements = computeDofsElements(obj,dofsEdges)
            mesh = obj.mesh;
            if obj.polynomialOrder == 1
                dofsElements = [];
            else
                ord = obj.polynomialOrder-2;
                ndofsElements = 0;
                for i = 0:ord
                    ndofsElements = ndofsElements + i;
                end
                dofsElements = zeros(obj.mesh.nelem,ndofsElements);
                for i = 1:obj.mesh.nelem
                    dofsElements(i,:) = (i-1)*ndofsElements+1:i*ndofsElements;
                end
                dofsElements = dofsElements + max(max(dofsEdges));
            end
        end
        
        function computeCoords(obj)
            nelem = size(obj.dofs,1);
            coord = obj.mesh.coord;
            
            x = zeros(obj.ndofs,1);
            y = zeros(obj.ndofs,1);
            
            if obj.interpolation.lagrangeElement.polynomialOrder~=1
                for ielem = 1:nelem
                    coor = obj.computeNodesElement(coord(obj.dofs(ielem,1:3),:));
                    x(obj.dofs(ielem,:)) = coor(:,1);
                    y(obj.dofs(ielem,:)) = coor(:,2);
                end
            
                obj.coords = [x,y];
            else
                obj.coords = obj.mesh.coord;
            end
        end
        
        function coor = computeNodesElement(obj,coords)
            base = obj.interpolation.lagrangeElement.nodes;          
            c = base(1:3,:);
            M = (coords-coords(1,:))'/c';
            N = coords(1,:)';
            coor = (M*base'+N)';
        end
        
    end
end