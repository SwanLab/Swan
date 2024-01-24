classdef P1Function < FeFunction

    properties (GetAccess = public, SetAccess = private)
        nDofs
        nDofsElem
    end

    properties (Access = private)
        interpolation
    end

    properties (Access = private)
        meshCoarse
        meshFine
    end

    methods (Access = public)

        function obj = P1Function(cParams)
            obj.init(cParams);
            obj.createInterpolation();
            obj.computeNDofs();
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

        function p1sub = evaluateSubdomain(obj, domain)
            mesh_sub = obj.mesh.getSubmesh(domain);
            dofs = domain(obj.mesh.coord);
            s.fValues = obj.fValues(dofs, :);
            s.mesh    = mesh_sub;
            p1sub = P1Function(s);
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
            nElem = size(grad.fValues, 3);
            nGauss = size(grad.fValues, 2);

            gradReshp = reshape(grad.fValues, [nDims,nDimf,nGauss,nElem]);
            gradT = permute(gradReshp, [2 1 3 4]);
            symGrad = 0.5*(gradReshp + gradT);
            
            rshp = reshape(symGrad, [nDims*nDimf,nGauss,nElem]);
            s.fValues = rshp;
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
            fFine = P1Function(s);
        end

        function dofConnec = computeDofConnectivity(obj)
            conne  = obj.mesh.connec;
            nDimf  = obj.ndimf;
            nNode  = size(conne, 2);
            nDofsE = nNode*nDimf;
            dofsElem  = zeros(nDofsE,size(conne,1));
            for iNode = 1:nNode
                for iUnkn = 1:nDimf
                    idofElem   = nDimf*(iNode - 1) + iUnkn;
                    globalNode = conne(:,iNode);
                    idofGlobal = nDimf*(globalNode - 1) + iUnkn;
                    dofsElem(idofElem,:) = idofGlobal;
                end
            end
            dofConnec = dofsElem;
        end

        function plot(obj) % 2D domains only

            switch obj.mesh.type
                case {'TRIANGLE','QUAD'}
                    x = obj.mesh.coord(:,1);
                    y = obj.mesh.coord(:,2);
                    figure()
                    for idim = 1:obj.ndimf
                        subplot(1,obj.ndimf,idim);
                        z = obj.fValues(:,idim);
                        a = trisurf(obj.mesh.connec,x,y,z);
                        view(0,90)
                        %             colorbar
                        shading interp
                        a.EdgeColor = [0 0 0];
                        title(['dim = ', num2str(idim)]);
                        colorbar
                    end
                case 'LINE'
                    x = obj.mesh.coord(:,1);
                    y = obj.fValues;
                    figure()
                    plot(x,y)
            end
        end

        function plotArrowVector(obj)
     %       figure()
            a = obj.fValues;
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            ax = a(:,1);
            ay = a(:,2);
            q = quiver(x,y,ax,ay);
            q.ShowArrowHead = 'off';
        end

        function print(obj, filename, software)
            if nargin == 2; software = 'GiD'; end
            s.mesh = obj.mesh;
            s.fun = {obj};
            s.type = software;
            s.filename = filename;
            p = FunctionPrinter.create(s);
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

        function v = computeL2norm(obj)
            s.type     = 'ScalarProduct';
            s.quadType = 'QUADRATICMASS';
            s.mesh     = obj.mesh;
            int = Integrator.create(s);
            ff  = int.compute(obj,obj);
            v   = sqrt(ff);
        end


    end

    methods (Access = public, Static)

        function p1 = create(mesh, ndimf)
            s.fValues = zeros(mesh.nnodes, ndimf);
            s.mesh    = mesh;
            p1 = P1Function(s);
        end

        function fS = times(f1,f2)
            fS = f1.fValues.*f2.fValues;
            s.fValues = fS;
            s.mesh    = f1.mesh;
            fS = P1Function(s);
        end

        function fS = sum(f1,f2)
            fS = f1.fValues+f2.fValues;
            s.fValues = fS;
            s.mesh    = f1.mesh;
            fS = P1Function(s);
        end
        
        
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
            obj.fValues = cParams.fValues;
            obj.ndimf   = size(cParams.fValues,2);
            obj.order   = 'LINEAR';
            obj.meshCoarse = cParams.mesh;
        end

        function createInterpolation(obj)
            m.type = obj.mesh.type;
            obj.interpolation = Interpolation.create(m,'LINEAR');
        end

        function computeNDofs(obj)
            obj.nDofsElem = obj.ndimf*obj.interpolation.nnode;
            obj.nDofs = obj.ndimf * size(obj.fValues, 1);
        end

        % Printing
        function fM = getFormattedFValues(obj)
            fM = obj.fValues;
        end

        function f = computeFunctionInEdges(obj,m,fNodes)
            s.edgeMesh = m.computeEdgeMesh();
            s.fNodes   = fNodes;
            eF         = EdgeFunctionInterpolator(s);
            f = eF.compute();
        end

    end

end
