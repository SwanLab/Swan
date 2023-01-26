classdef P0Function < FeFunction
    
    properties (Access = public)
    end

    properties (Access = private)
        connec
        type
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

        function fDfun = computeP1DiscontinuousFunction(obj)
            fD = obj.createDiscontinuousFunction();
            s.fValues = fD;
            s.connec  = obj.connec;
            s.type    = obj.type;
            fDfun = P1DiscontinuousFunction(s);            
        end

        function plot(obj, m)
            p1DiscFun = obj.computeP1DiscontinuousFunction(); % replace with a projector
            p1DiscFun.plot(m);
        end

        function print(obj, s)
%             s.mesh
            s.fun = {obj};
            p = FunctionPrinter(s);
            p.print();
        end

        function [res, pformat] = getDataToPrint(obj)
            pformat = ['%s ',repmat('%12.5d ',1,obj.ndimf),'\n'];
            elemColum = obj.computeElementStringColum();
            valColums = obj.computeTensorValueColums();
            res = [elemColum,valColums]';
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fValues = cParams.fValues;
            obj.connec  = cParams.connec; % Needed for discontinuous
            obj.type    = cParams.type;
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
            nnodeElem = size(obj.connec,2);
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
        function c = computeElementStringColum(obj)
            q = Quadrature.set(obj.type);
            q.computeQuadrature('LINEAR');
            nElem = size(obj.connec, 1);
            nGaus = q.ngaus;
            allElem(:,1) = 1:nElem;
            colWidth = size(num2str(nElem),2);
            strInCol = repmat(' ',nElem*nGaus,colWidth);
            numIndex = 1:nGaus:nElem*nGaus;
            strInCol(numIndex,:) = num2str(allElem);
            c = cellstr(strInCol);
        end
        
        function fM = computeTensorValueColums(obj)
%             fV = obj.fValues;
            q = Quadrature.set(obj.type);
            q.computeQuadrature('LINEAR');
            fV = obj.evaluate(q.posgp);
            nGaus   = q.ngaus;
            nComp   = obj.ndimf;
            nElem   = size(obj.connec, 1);
            fM  = zeros(nGaus*nElem,nComp);
            for iStre = 1:nComp
                for iGaus = 1:nGaus
                    rows = linspace(iGaus,(nElem - 1)*nGaus + iGaus,nElem);
                    fM(rows,iStre) = fV(iStre,iGaus,:);
                end
            end
            fM = num2cell(fM);
        end

    end

end