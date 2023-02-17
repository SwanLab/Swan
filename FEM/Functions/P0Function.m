classdef P0Function < FeFunction
    
    properties (Access = public)
    end
    
    methods (Access = public)
        
        function obj = P0Function(cParams)
            obj.init(cParams);
            obj.createFvaluesByElem();
        end

        function fxV = evaluate(obj, xV)
            % Its a p0 function, so no true need to interpolate -- the
            % value is constant
            nGaus = size(xV,2);
            nFlds = size(obj.fValues,1);
            nElem = size(obj.fValues,3);
            fxV = zeros(nFlds,nGaus,nElem);
            for iGaus = 1:nGaus
                fxV(:,iGaus,:) = squeeze(obj.fValues);
            end
        end

        function fD = computeP1DiscontinuousFunction(obj)
            s.mesh    = obj.mesh;
            s.connec  = obj.mesh.connec;
            p = Projector_toP1Discontinuous(s);
            fD = p.project(obj);
        end

        function plot(obj)
            p1DiscFun = obj.computeP1DiscontinuousFunction(); % replace with a projector
            p1DiscFun.plot();
        end

        function print(obj, s)
            s.mesh = obj.mesh;
            s.fun = {obj};
            p = FunctionPrinter(s);
            p.print();
        end

        function [res, pformat] = getDataToPrint(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            nElem = size(obj.mesh.connec, 1);
            nGaus = q.ngaus;

            s.nDimf   = obj.ndimf;
            s.nData   = nElem*nGaus;
            s.nGroup  = nElem;
            s.fValues = obj.getFormattedFValues();
            fps = FunctionPrintingSettings(s);
            [res, pformat] = fps.getDataToPrint();
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fValues = cParams.fValues;
            obj.mesh    = cParams.mesh;
            obj.ndimf   = size(cParams.fValues,2);
        end

        function createFvaluesByElem(obj)
            f = obj.fValues;
            nElem = size(f,1);
            nDime = size(f,2);
            obj.fValues = reshape(f',[nDime, 1, nElem]);
        end

        function fD = createDiscontinuousFunction(obj)
            dim = 1;
            ndim  = size(obj.fValues, 1);
            nnodeElem = size(obj.mesh.connec,2);
            fEl = squeeze(obj.fValues(dim,:,:));
            fRepeated = zeros(ndim, size(fEl,1), nnodeElem);
            for idim = 1:ndim
                fEl = squeeze(obj.fValues(idim,:,:));
                for iNode = 1:nnodeElem
                    fRepeated(idim, :,iNode) = fEl;
                end
            end
            fD = permute(fRepeated, [1 3 2]);
        end

        % Printing
        function fM = getFormattedFValues(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            fV = obj.evaluate(q.posgp);
            nGaus   = q.ngaus;
            nComp   = obj.ndimf;
            nElem   = size(obj.mesh.connec, 1);
            fM  = zeros(nGaus*nElem,nComp);
            for iStre = 1:nComp
                for iGaus = 1:nGaus
                    rows = linspace(iGaus,(nElem - 1)*nGaus + iGaus,nElem);
                    fM(rows,iStre) = fV(iStre,iGaus,:);
                end
            end
        end

    end

end