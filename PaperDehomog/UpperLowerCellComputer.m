classdef UpperLowerCellComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
        isUpperCell
        isLowerCell
    end
    
    properties (Access = private)
        isNotCoherent
        isCoherent
        referenceCells
        isVertexInCell
    end
    
    properties (Access = private)
        mesh
        areCoherent
        pathVertexes
    end
    
    methods (Access = public)
        
        function obj = UpperLowerCellComputer(cParams)
            obj.init(cParams)
        end
        
        function compute(obj)
            obj.computeReferenceCells();
            for ivertex = 1:length(obj.pathVertexes)
                obj.computeIsVertexInCell(ivertex);
                obj.isCellAnUpperCell(ivertex);
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.pathVertexes = cParams.pathVertexes;
            obj.areCoherent  = cParams.areCoherent;
            obj.isUpperCell  = obj.createEmptyP0Function();
            obj.isLowerCell  = obj.createEmptyP0Function();
        end
        
        function computeReferenceCells(obj)
            s.mesh         = obj.mesh;
            s.pathVertexes = obj.pathVertexes;
            r = ReferenceCellOfPathComputer(s);
            rC = r.compute();
            obj.referenceCells = rC;
        end
        
        function computeIsVertexInCell(obj,ivertex)
            vI = obj.pathVertexes(ivertex);
            vertexInCell  = obj.mesh.connec;
            itIs = (vertexInCell == vI);
            obj.isVertexInCell = itIs;
        end
        
        function isCellAnUpperCell(obj,ivertex)
            rCell  = obj.referenceCells(ivertex);
            obj.computeCoherentAndNotCoherentVertex();
            [itHas,itHasNot] = obj.hasCellSameCoherentOrientationAsReferenceCell(rCell);
            isU = squeeze(obj.isUpperCell.fValues);
            isL = squeeze(obj.isLowerCell.fValues);
%             if rCell == obj.singularElement
%                 isU(rCell) = true;
%             end
            isU(itHas)    = isU(rCell);
            isU(itHasNot) = ~isU(rCell);

            isL(itHas)    = ~isU(rCell);
            isL(itHasNot) = isU(rCell);
            obj.isUpperCell.setFValues(isU);
            obj.isLowerCell.setFValues(isL);
        end
        
        function computeCoherentAndNotCoherentVertex(obj)
            areC  = squeeze(obj.areCoherent.fValues(1,:,:))';
            isCoh = obj.restrictToCell(areC);
            isNot = obj.restrictToCell(~areC);
            obj.isCoherent    = isCoh;
            obj.isNotCoherent = isNot;
        end

        function [itHas,itHasNot] = hasCellSameCoherentOrientationAsReferenceCell(obj,rCell)
            if obj.isCoherent(rCell)
                itHas    = obj.isCoherent;
                itHasNot = obj.isNotCoherent;
            else
                itHas    = obj.isNotCoherent;
                itHasNot = obj.isCoherent;
            end
        end
        
        function f = restrictToCell(obj,f)
            fV = obj.restrictToVertex(f);
            f  =  any(fV,2);
        end
        
        function fV = restrictToVertex(obj,f)
            itIs = obj.isVertexInCell;
            fV = f & itIs;
        end
        
        function f0 = createEmptyP0Function(obj)
            s.mesh    = obj.mesh;
            s.order   = 'P0';
            s.fValues = false(obj.mesh.nelem,1);
            f0 = LagrangianFunction(s);
        end

    end
    
end