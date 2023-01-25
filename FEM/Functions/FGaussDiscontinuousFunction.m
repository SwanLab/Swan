 classdef FGaussDiscontinuousFunction < handle
    
    properties (Constant, Access = public)
        fType = 'GAUSSPOINTS'
    end

    properties (Access = public)
        ndimf
        fValues
        quadrature
    end

    properties (Access = private)
    end
    
    properties (Access = private)
    end
    
    methods (Access = public)
        
        function obj = FGaussDiscontinuousFunction(cParams)
            obj.init(cParams)
        end

        function fxV = evaluate(obj, xV)
            assert(isequal(xV, obj.quadrature.posgp), 'Gauss points do not match')
            fxV = obj.fValues;
        end

        function plot(obj, mesh)
            pp1.mesh   = mesh;
            pp1.connec = mesh.connec;
            projP1 = Projector_toP1(pp1);
            p1fg = projP1.project(obj);
            p1fg.plot(mesh);
        end

        function print(obj, s)
%             s.mesh
            s.fun = obj;
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
            obj.fValues    = cParams.fValues;
            obj.quadrature = cParams.quadrature;
            obj.ndimf      = size(cParams.fValues,1);
        end

        % Printing
        function c = computeElementStringColum(obj)
            nElem = size(obj.fValues, 3);
            nGaus = obj.quadrature.ngaus;
            allElem(:,1) = 1:nElem;
            colWidth = size(num2str(nElem),2);
            strInCol = repmat(' ',nElem*nGaus,colWidth);
            numIndex = 1:nGaus:nElem*nGaus;
            strInCol(numIndex,:) = num2str(allElem);
            c = cellstr(strInCol);
        end
        
        function fM = computeTensorValueColums(obj)
            fV = obj.fValues;
            nGaus = obj.quadrature.ngaus;
            nComp = obj.ndimf;
            nElem = size(obj.fValues, 3);
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