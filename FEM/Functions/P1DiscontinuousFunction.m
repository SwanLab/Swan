classdef P1DiscontinuousFunction < FeFunction
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        interpolation
    end
    
    properties (Access = private)
    end
    
    methods (Access = public)
        
        function obj = P1DiscontinuousFunction(cParams)
            obj.init(cParams)
            obj.createInterpolation();
        end

        function fxV = evaluate(obj, xV)
            % Goal: develop this function
            func = obj.fValues;
            obj.interpolation.computeShapeDeriv(xV);
            shapes = obj.interpolation.shape;
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
        
        function fFine = refine(obj, m, mFine)
         %   mFineD = mFine.createDiscontinuousMesh();
            f = squeeze(obj.fValues);
            f = f(:);
            fEdges = obj.computeFunctionInEdges(m,f);
            fAll  = [f;fEdges];
            
            s.mesh    = mFine;
            s.fValues = fAll;
            p1fun = P1Function(s);
            fFine = p1fun.project('P1D');
        end

        function fV = getFvaluesAsVector(obj)
            ndims   = size(obj.fValues, 1);
            nelem   = size(obj.mesh.connec, 1);
            nnodeEl = size(obj.mesh.connec, 2);
            fV = reshape(obj.fValues, [ndims, nelem*nnodeEl])';
        end

        function plot(obj)
            fD = obj.getFvaluesAsVector();
            mD = obj.mesh.createDiscontinuousMesh();
            x = mD.coord(:,1);
            y = mD.coord(:,2);
            figure()
            for idim = 1:obj.ndimf
                subplot(1,obj.ndimf,idim);
                z = fD(:,idim);
                a = trisurf(mD.connec,x,y,z);
                view(0,90)
    %             colorbar
                shading interp
                a.EdgeColor = [0 0 0];
                title(['dim = ', num2str(idim)]);
            end
        end

        function print(obj, s)
            s.mesh = obj.mesh.createDiscontinuousMesh();
            s.fun = {obj};
            p = FunctionPrinter(s);
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

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fValues = cParams.fValues;
            obj.mesh    = cParams.mesh;
            obj.ndimf   = size(cParams.fValues,1);
        end

        function createInterpolation(obj)
            m.type = obj.mesh.type;
            obj.interpolation = Interpolation.create(m,'LINEAR');
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