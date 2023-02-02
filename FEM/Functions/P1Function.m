classdef P1Function < FeFunction

    properties (Access = public)
    end

    properties (Access = private)
        interpolation
    end

    properties (Access = private)
        type
        connec
    end

    methods (Access = public)

        function obj = P1Function(cParams)
            obj.init(cParams);
            obj.createInterpolation();
        end

        function fxV = evaluate(obj, xV)
            obj.interpolation.computeShapeDeriv(xV);
            shapes = obj.interpolation.shape;
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
        
        function dNdx  = computeCartesianDerivatives(obj,quad,mesh)
            nElem = size(obj.connec,1);
            nNode = obj.interpolation.nnode;
            nDime = obj.interpolation.ndime;
            nGaus = quad.ngaus;
            invJ  = mesh.computeInverseJacobian(quad,obj.interpolation);
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

        function gradFun = computeGradient(obj, quad, mesh)
            dNdx = obj.computeCartesianDerivatives(quad,mesh);
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
            s.quadrature = quad;
            gradFun = FGaussDiscontinuousFunction(s);
        end

        function symGradFun = computeSymmetricGradient(obj,quad,mesh)
            grad = obj.computeGradient(quad,mesh);
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

        function fFine = refine(obj,m,mFine)
            fNodes  = obj.fValues;
            fEdges  = obj.computeFunctionInEdges(m,fNodes);
            fAll    = [fNodes;fEdges];
            s.type    = mFine.type;
            s.connec  = mFine.connec;
            s.fValues = fAll;
            fFine = P1Function(s);
        end

        function plot(obj, m) % 2D domains only
            x = m.coord(:,1);
            y = m.coord(:,2);
            figure()
            for idim = 1:obj.ndimf
                subplot(1,obj.ndimf,idim);
                z = obj.fValues(:,idim);
                a = trisurf(m.connec,x,y,z);
                view(0,90)
    %             colorbar
                shading interp
                a.EdgeColor = [0 0 0];
                title(['dim = ', num2str(idim)]);
            end
        end

        function print(obj, s)
%             s.mesh
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

        function fD = createP1Discontinous(obj,m)
            s.mesh   = m;
            s.connec = obj.connec;
            p = ProjectorToP1discont(s);
            s.x = obj;
            s.origin = 'P1';
            fD = p.project(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.connec  = cParams.connec;
            obj.type    = cParams.type;
            obj.fValues = cParams.fValues;
            obj.ndimf   = size(cParams.fValues,2);
            obj.order   = 'LINEAR';
        end

        function createInterpolation(obj)
            m.type = obj.type;
            obj.interpolation = Interpolation.create(m,'LINEAR');
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