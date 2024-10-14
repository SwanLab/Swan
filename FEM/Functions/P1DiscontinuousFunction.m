classdef P1DiscontinuousFunction < FeFunction

    properties (Access = public)
    end

    properties (Access = private)
        fValuesDisc
        interpolation        
    end

    properties (Access = public)
        dofConnec
        dofCoord
    end

    methods (Access = public)

        function obj = P1DiscontinuousFunction(cParams)
            obj.init(cParams)
            obj.order = '1';
            obj.createInterpolation();
          %  obj.createDOFCoordConnec();
           obj.createValuesByElement();
        end

        function fxV = evaluate(obj, xV)
            func = obj.fValuesDisc;
            shapes = obj.interpolation.computeShapeFunctions(xV);
            nNode  = size(shapes,1);
            nGaus  = size(shapes,2);
            nF     = size(func,1);
            nElem  = size(func,3);
            fxV = zeros(nF,nGaus,nElem);
            for kNode = 1:nNode
                shapeKJ = shapes(kNode,:,:);
                fKJ     = func(:,kNode,:);
                f = bsxfun(@times,shapeKJ,fKJ);
                fxV = fxV + f;
            end
        end

        function f = copy(obj)
            f = obj.create(obj.mesh, obj.dofConnec, obj.dofCoord,obj.ndimf);
            f.fValues = obj.fValues;
        end        

        function ord = getOrderNum(obj)
            ord = str2double(obj.order(end));
        end

        function c = getDofConnec(obj)
            c = obj.dofConnec;
        end

        function c = getDofCoord(obj)
            c = obj.dofCoord;
        end

        function fD = getFvaluesDisc(obj)
            obj.createValuesByElement();            
            fD = obj.fValuesDisc;
        end

        function N = computeShapeFunctions(obj, xV)
            N = obj.interpolation.computeShapeFunctions(xV);
        end

        function dN = computeShapeDerivatives(obj,xV)
            dN = obj.interpolation.computeShapeDerivatives(xV);
        end

        function dNdx  = evaluateCartesianDerivatives(obj,xV)
            nElem = size(obj.dofConnec,1);
            nNodeE = obj.interpolation.nnode;
            nDimE = obj.interpolation.ndime;
            nDimG = obj.mesh.ndim;
            nPoints = size(xV, 2);
            invJ  = obj.mesh.computeInverseJacobian(xV);
            deriv = obj.computeShapeDerivatives(xV);
            dShapes  = zeros(nDimG,nNodeE,nPoints,nElem);
            for iDimG = 1:nDimG
                for kNodeE = 1:nNodeE
                    for jDimE = 1:nDimE
                        invJ_IJ   = invJ(iDimG,jDimE,:,:);
                        dShapes_JK = deriv(jDimE,kNodeE,:);
                        dShapes_KI   = pagemtimes(invJ_IJ,dShapes_JK);
                        dShapes(iDimG,kNodeE,:,:) = dShapes(iDimG,kNodeE,:,:) + dShapes_KI;
                    end
                end
            end
            dNdx = dShapes;
        end

        function gradFun = evaluateGradient(obj, xV)
            dNdx = obj.evaluateCartesianDerivatives(xV);
            nDimf = obj.ndimf;
            nDims = size(dNdx, 1); % derivX, derivY (mesh-related?)
            nNode = size(dNdx, 2);
            nElem = size(dNdx, 3);
            nGaus = size(dNdx, 4);


            cV = squeeze(obj.fValues);
            grad = zeros(nDims,nDimf, nElem, nGaus);
            for iGaus = 1:nGaus
                dNdx_g = dNdx(:,:,:,iGaus);
                for iDims = 1:nDims
                    for iNode = 1:nNode
                        dNdx_i = squeeze(dNdx_g(iDims, iNode,:));
                        %dNi = squeeze(dNdx(idim,iNode,:,igaus));
                        f = cV(iNode,:)';
                        gradV = dNdx_i.*f;
                        fG(:,1) = squeeze(grad(iDims,:,:,iGaus));
                        grad(iDims,:,:,iGaus) = fG + gradV;
                    end
                end
            end
            %  gradt(1,:,:,:) = grad;
            fVR = reshape(grad, [nDims*nDimf,nElem, nGaus]);
            s.fValues = permute(fVR, [1 3 2]);
            s.mesh    = obj.mesh;
            s.quadrature = xV;
            gradFun = FGaussDiscontinuousFunction(s);
        end

        function fFine = refine(obj,fD, mFine)
            P1Dref = P1Refiner(fD);
            fFine  = P1Dref.compute();

            f = obj.fValuesDisc;
            for iDim = 1:obj.ndimf
                fI = f(iDim,:,:);
                fI = fI(:);
                fEdges = obj.computeFunctionInEdges(obj.mesh,fI);
                fAll(:,iDim)  = [fI;fEdges];
            end
            s.mesh    = mFine;
            s.fValues = fAll;
            s.order   = 'P1';
            p1fun = LagrangianFunction(s);
            fFine = p1fun.project('P1D');
        end

        % function dofConnec = computeDofConnectivity(obj)
        %     conne  = obj.dofConnec;
        %     nDimf  = obj.ndimf;
        %     nNodeE = size(conne, 2);
        %     nDofsE = nNodeE*nDimf;
        %     dofsElem  = zeros(nDofsE,size(conne,1));
        %     for iNode = 1:nNodeE
        %         for iUnkn = 1:nDimf
        %             idofElem   = nDimf*(iNode - 1) + iUnkn;
        %             globalNode = conne(:,iNode);
        %             idofGlobal = nDimf*(globalNode - 1) + iUnkn;
        %             dofsElem(idofElem,:) = idofGlobal;
        %         end
        %     end
        %     dofConnec = dofsElem;
        % end

        function fR = getFvaluesAsVector(obj)
            f  = obj.fValues;
            fR = obj.reshapeAsVector(f);
        end

        function fR = reshapeAsVector(obj,fValues)
            ndims   = size(fValues, 1);
            nelem   = size(obj.mesh.connec, 1);
            nnodeEl = size(obj.mesh.connec, 2);
            fR = reshape(fValues, [ndims, nelem*nnodeEl])';
        end

        function plot(obj)
            for iDim = 1:obj.ndimf
                connecf{iDim} = obj.getDofConnecByVector();
                coordf{iDim}  = obj.getDofCoordByVector(iDim);
            end
            s.connec = connecf;
            s.coord  = coordf;



        
            s.fValues = obj.fValues;
            s.ndimf   = obj.ndimf;
            lP = LagrangianPlotter(s);
            lP.plot();
        end

        function node = getDofConnecByVector(obj)
            nNode = obj.interpolation.nnode;
            for iNode = 1:nNode
                iDof   = (iNode-1)*obj.ndimf+1;
                node(:,iNode) = (obj.dofConnec(:,iDof)-1)/obj.ndimf+1;
            end
        end

        function cT = getDofCoordByVector(obj,dimf)
            for iDim = 1:obj.mesh.ndim
                coordN = obj.dofCoord(:,iDim);
                cResh  = reshape(coordN',obj.ndimf,[]);
                cT(:,iDim) = cResh(dimf,:);
            end
        end

        function plotLine(obj)
            y = obj.mesh.coord(:,2);
            nelem = size(obj.fValues,3);
            ndim  = size(obj.fValues,1);
            for idim=1:ndim
                yPlot = []; xPlot = [];
                for e=1:nelem
                    yPlot = [yPlot y(e,1) y(e+1,1)];
                    xPlot = [xPlot obj.fValues(idim,:,e)];

                end
                plot(xPlot,yPlot)
                hold on
            end
        end

        function plotContour(obj)
            fD = obj.getFvaluesAsVector();
            xy = obj.reshapeAsVector(obj.dofCoord);            
            x = xy(:,1);
            y = xy(:,2);
            figure()
            for idim = 1:obj.ndimf
                subplot(1,obj.ndimf,idim);
                z = fD(:,idim);
                [~,a] = tricontour(obj.dofConnec,x,y,z,30);
                set(a,'LineWidth',5);
                view(0,90)
                colorbar
                title(['dim = ', num2str(idim)]);
            end
        end

        function print(obj, filename, software)
            if nargin == 2; software = 'Paraview'; end
            s.mesh = obj.mesh;
            s.fun = {obj};
            s.type = software;
            s.filename = filename;
            p = FunctionPrinter.create(s);
            p.print();
        end

        function [res, pformat] = getDataToPrint(obj)
            nElem = size(obj.fValues, 3);
            nNodE = size(obj.fValues, 2);
            nNods = nElem*nNodE;
            s.nDimf   = obj.ndimf;
            s.nData   = nNods;
            s.nGroup  = nNods;
            s.fValues = obj.getFormattedFValues();
            fps = FunctionPrintingSettings(s);
            [res, pformat] = fps.getDataToPrint();
        end

        function createDOFConnec(obj)
            nNodes = obj.mesh.nnodeElem*obj.mesh.nelem;
            nodes  = 1:nNodes;
            c      = reshape(nodes,obj.mesh.nnodeElem,obj.mesh.nelem)';
            obj.dofConnec = c;
        end

        function createDOFCoord(obj)
            conn   = obj.mesh.connec;
            coordC = obj.mesh.coord;
            nNode = size(conn,2);
            nDime = size(coordC,2);
            nodes = reshape(conn',1,[]);
            fe    = coordC(nodes,:)';
            coorD = reshape(fe,nDime,nNode,[]);
            obj.dofCoord = coorD;
        end

        function createDOFCoordConnec(obj)
            obj.createDOFCoord();
            obj.createDOFConnec()
        end

        function ord = orderTextual(obj)
            ord = 'LINEAR';
        end

    end


    methods (Access = public, Static)

        function fS = sum(f1,f2)
            fS = f1.fValues + f2.fValues;
            s.fValues = fS;
            s.mesh    = f1.mesh;
            fS = P1DiscontinuousFunction(s);
        end

        function p1d = create(mesh, dofConnec, dofCoord,ndimf)
            s.mesh      = mesh;
            s.dofConnec = dofConnec;
            s.dofCoord  = dofCoord;  
            s.fValues = zeros(mesh.nnodeElem*mesh.nelem,ndimf);
            p1d = P1DiscontinuousFunction(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh      = cParams.mesh;            
            obj.fValues   = cParams.fValues;
            obj.dofConnec = cParams.dofConnec;
            obj.dofCoord  = cParams.dofCoord;
            obj.ndimf     = size(cParams.fValues,2);
            obj.order     = 'LINEAR';
        end

        function createValuesByElement(obj)
            nDimF  = size(obj.fValues,2);            
            fVals = reshape(obj.fValues',nDimF,[],obj.mesh.nelem);
            obj.fValuesDisc = fVals;
        end

        function createInterpolation(obj)
            type = obj.mesh.type;
            obj.interpolation = Interpolation.create(type,'LINEAR');
        end

        % Printing
        function fM = getFormattedFValues(obj)
            nComp = obj.ndimf;
            nNodE = size(obj.fValues, 2);
            nElem = size(obj.fValues, 3);
            nVals = nNodE*nElem;
            fM = reshape(obj.fValues, [nComp, nVals])';
        end

        function f = computeFunctionInEdges(obj,m,fNodes)
            s.edgeMesh = m.computeEdgeMesh();
            s.fNodes   = fNodes;
            eF         = EdgeFunctionInterpolator(s);
            f = eF.compute();
        end

    end


end
