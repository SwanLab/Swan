classdef SubCellsCasesComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        caseInfo
    end
    
    properties (Access = private)
        isNodeInterior
        intergerCodeCases
        integerCases
    end
    
    properties (Access = private)
        connec
        levelSet
    end
    
    methods (Access = public)
        
        function obj = SubCellsCasesComputer(cParams)
            obj.init(cParams);
            obj.computeCutCase();
            obj.computeIntegerCases();
        end
        %
        function c = computeSubCellsTriangle(obj)
            nElem  = size(obj.isNodeInterior,1);
            nCases = size(obj.intergerCodeCases,1);
            subCell = false(nElem,nCases);
            nSubCells = 3;
            isInt = false(nSubCells,nElem);
            for icase = 1:nCases
                isCaseA = obj.computeIntegerCase(icase,1);
                isCaseB = obj.computeIntegerCase(icase,2);
                isCase = or(isCaseA,isCaseB);
                subCell(:,icase) = isCase;
                isInt(2:end,isCaseA) = true;
                isInt(1,isCaseB) = true;
            end 
  %          c.cellsCase          = true(size(obj.isNodeInterior,1));
            c.subCellCases       = subCell;
            c.isSubCellsInterior  = isInt;
        end
        
        function c = computeSubCellsThetaedraThreeVsOne(obj)
            nElem  = size(obj.isNodeInterior,1);
            nCases = size(obj.intergerCodeCases,1);
            subCell = false(nElem,nCases);
            nSubCells = 4;
            isInt = false(nSubCells,nElem);
            for icase = 1:nCases
                isCaseA = obj.computeIntegerCase(icase,1);
                isCaseB = obj.computeIntegerCase(icase,2);
                isCase = or(isCaseA,isCaseB);
                subCell(:,icase) = isCase;
                isInt(2:end,isCaseA) = true;
                isInt(1,isCaseB) = true;
            end
            c.cellsCase          = obj.isThreeVsOne();
            c.subCellCases       = subCell;
            c.isSubCellsInterior  = isInt;
        end
        
        function c = computeSubCellsThetaedraTwoVsTwo(obj)
            nElem  = size(obj.isNodeInterior,1);
            nCases = size(obj.intergerCodeCases,1);
            subCell = false(nElem,nCases);
            nSubCells = 6;
            isInt = false(nSubCells,nElem);
            for icase = 1:nCases
                isCaseA = obj.computeIntegerCase(icase,1);
                isCaseB = obj.computeIntegerCase(icase,2);
                isCase = or(isCaseA,isCaseB);
                subCell(:,icase) = isCase;
                isInt(1:3,isCaseA) = true;
                isInt(4:6,isCaseB) = true;
            end
            c.cellsCase          = obj.isTwoVsTwo();
            c.subCellCases       = subCell;
            c.isSubCellsInterior  = isInt;
        end
        
        
        function compute(obj)
                c = cell(3,1);
            switch size(obj.isNodeInterior,2)
                case 3
                    obj.intergerCodeCases = [6 1;5 2;3 4];
                    c{1} = obj.computeSubCellsTriangle();
                case 4                    
                  %  switch mode(sum(obj.isNodeInterior,2))
                       % case {1,3}
                            obj.intergerCodeCases = [14 1;13 2;11 4;7 8];
                            c{2} = obj.computeSubCellsThetaedraThreeVsOne();
                       % case 2
                            obj.intergerCodeCases = [9 6;3 12;5 10];
                            c{3} = obj.computeSubCellsThetaedraTwoVsTwo();
                 %   end
            end
            obj.caseInfo = c;
        end
        
        function itIs = isThreeVsOne(obj)
            nInteriorNode = sum(obj.isNodeInterior,2);
            isCase1 = nInteriorNode == 1;
            isCase3 = nInteriorNode == 3;
            itIs = or(isCase1,isCase3);
        end
        
        function itIs = isTwoVsTwo(obj)
            nInteriorNode = sum(obj.isNodeInterior,2);
            itIs = nInteriorNode == 2;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.connec  = cParams.connec;
            obj.levelSet = cParams.levelSet;
        end
        
        
        function computeCutCase(obj)
            nodes = obj.connec;
            ls = zeros(size(nodes));
            for iNode = 1:size(nodes,2)
                ls(:,iNode) = obj.levelSet(nodes(:,iNode));
            end
            obj.isNodeInterior = 1 - heaviside(ls);
        end
        
        function d = computeIntegerCases(obj)
            nodes = obj.isNodeInterior;
            nnode = size(nodes,2);
            nodePos = (1:nnode) - 1;
            pow2vector = 2.^(nodePos);
            d = pow2vector*nodes';
            obj.integerCases = d;
        end
        
        function isCase = computeIntegerCase(obj,icase,subCase)
            integerCase = obj.integerCases;
            isCase = integerCase == obj.intergerCodeCases(icase,subCase);
        end
        
        
        
    end
    
end