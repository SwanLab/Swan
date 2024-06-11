classdef P2Function < FeFunction

    properties (Access = public)
    end

    properties (Access = private)
        interpolation
        coord
        connec
    end

    methods (Access = public)

        function obj = P2Function(cParams)
            obj.init(cParams);
            obj.createInterpolation();
            obj.createDOFCoordConnec();
        end

        function fxV = evaluate(obj, xV)
            shapes = obj.interpolation.computeShapeFunctions(xV);
            nNode  = size(shapes,1);
            nGaus  = size(shapes,2);
            nF     = size(obj.fValues,2);
            nElem  = size(obj.connec,1);
            fxV = zeros(nF,nGaus,nElem);
            for iGaus = 1:nGaus
                for iNode = 1:nNode
                    node = obj.connec(:,iNode);
                    Ni = shapes(iNode,iGaus);
                    fi = obj.fValues(node,:);
                    f(:,1,:) = Ni*fi';
                    fxV(:,iGaus,:) = fxV(:,iGaus,:) + f;
                end
            end

        end

        function N = computeShapeFunctions(obj, xV)
            N = obj.interpolation.computeShapeFunctions(xV);
        end
        
        function dNdx  = evaluateCartesianDerivatives(obj,xV)
            nElem = size(obj.connec,1);
            nNode = obj.interpolation.nnode;
            nDime = obj.interpolation.ndime;
            nGaus = size(xV,2);
            invJ  = obj.mesh.computeInverseJacobian(xV);
            deriv = obj.interpolation.computeShapeDerivatives(xV);
            dShapeDx  = zeros(nDime,nNode,nElem,nGaus);
            for igaus = 1:nGaus
                dShapes = deriv(:,:,igaus);
                for jDime = 1:nDime
                    invJ_JI   = invJ(:,jDime,:,igaus);
                    dShape_KJ = dShapes(jDime,:);
                    dSDx_KI   = bsxfun(@times, invJ_JI,dShape_KJ);
                    dShapeDx(:,:,:,igaus) = dShapeDx(:,:,:,igaus) + dSDx_KI;
                end
            end
            dNdx = dShapeDx;
        end

        function gradFun = computeGradient(obj, xV)
            dNdx = obj.evaluateCartesianDerivatives(xV);
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
                        nodes = obj.connec(:,iNode);
                        f = obj.fValues(nodes,:);
                        p = (dNdx_i.*f)';
                        pp(1,:,:) = p;
                        grad(iDims,:,:,iGaus) = grad(iDims,:,:,iGaus) + pp;
                    end
                end
            end
            fVR = reshape(grad, [nDims*nDimf,nElem, nGaus]);
            s.fValues = permute(fVR, [1 3 2]);
%             s.ndimf      = nDimf;
            s.quadrature = xV;
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
            symGradFun = FGaussDiscontinuousFunction(s);
        end

        function plot(obj) % 2D domains only
            s.mesh          = obj.mesh;
            s.interpolation = obj.interpolation;
            c = ConnecCoordFromInterpAndMesh(s);
            c.compute();
            coord = c.coord;
            connec = obj.connec(:, [1 4 2 5 3 6]);
            x = coord(:,1);
            y = coord(:,2);
            figure()
            for idim = 1:obj.ndimf
                subplot(1,obj.ndimf,idim);
                z = obj.fValues(:,idim);
                a = trisurf(connec,x,y,z);
                view(0,90)
    %             colorbar
                shading interp
                a.EdgeColor = [0 0 0];
                title(['dim = ', num2str(idim)]);
            end
        end

        function dofConnec = computeDofConnectivity(obj)
            conne  = obj.connec;
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

        function dof = getDofsFromCondition(obj, condition)
            nodes = condition(obj.coord);
            iNode = find(nodes==1);
            dofElem = repmat(1:obj.ndimf, [length(iNode) 1]);
            dofMat = obj.ndimf*(iNode - 1) + dofElem;
            dof = sort(dofMat(:));
        end

        function nod = getNodesFromCondition(obj, condition)
            nodes = condition(obj.coord);
            iNode = find(nodes==1);
%             dofElem = repmat(1:obj.ndimf, [length(iNode) 1]);
%             dofMat = obj.ndimf*(iNode - 1) + dofElem;
            nod = sort(iNode(:));
        end

    end

    methods (Access = public, Static)

        function p2 = create(mesh, ndimf)
            s.fValues = zeros(mesh.nnodes, ndimf); % wrong
            s.mesh    = mesh;
            p2 = P2Function(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.fValues = cParams.fValues;
            obj.ndimf   = size(cParams.fValues,2);
            obj.order   = 'ORDER4';
        end

        function createInterpolation(obj)
            type = obj.mesh.type;
            obj.interpolation = Interpolation.create(type,'QUADRATIC');
        end

        function createDOFCoordConnec(obj)
            s.mesh          = obj.mesh;
            s.interpolation = obj.interpolation;
            c = ConnecCoordFromInterpAndMesh(s);
            c.compute();
            obj.coord  = c.coord;
            obj.connec = c.connec;
            nDimf = size(obj.fValues,2);
            
            if all(obj.fValues == 0)
                obj.fValues = zeros(size(c.coord,1),nDimf);
            end
            
            if isequal(size(obj.mesh.coord,1), size(obj.fValues,1))
                obj.fValues = zeros(size(obj.coord,1),nDimf);
            end
        end

    end

end