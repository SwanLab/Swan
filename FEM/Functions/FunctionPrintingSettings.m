classdef FunctionPrintingSettings < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        nDimf
        nData
        nGroup
        fValues
    end
    
    methods (Access = public)

        function obj = FunctionPrintingSettings(cParams)
            obj.init(cParams);
        end

        function [res, pformat] = getDataToPrint(obj)
            nDimf = obj.nDimf;
            pformat = ['%s ',repmat('%12.5d ', 1, nDimf),'\n'];
            elemColum = obj.computeElementStringColum();
            valColums = obj.computeTensorValueColums();
            res = [elemColum,valColums]';
        end
        
    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.nData  = cParams.nData;
            obj.nGroup = cParams.nGroup;
            obj.nDimf  = cParams.nDimf;
            obj.fValues = cParams.fValues; % formatted for printing
        end

        function c = computeElementStringColum(obj)
            nDat = obj.nData;
            nGrp = obj.nGroup;
            nStp = nDat/nGrp;

            allElem(:,1) = 1:nGrp;
            colWidth = size(num2str(nGrp),2);
            strInCol = repmat(' ', nDat, colWidth);
            numIndex = 1:nStp:nDat;
            strInCol(numIndex,:) = num2str(allElem);
            c = cellstr(strInCol);
        end
        
        function fM = computeTensorValueColums(obj)
            fV = obj.fValues;
            fM = num2cell(fV);
        end

    end
    
end

