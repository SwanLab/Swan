classdef CorrectorComputer < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        correctorValues
        referenceCells
        isCoherent
        isNotCoherent
        itHasSameCoherence
        itHasNotSameCoherence
        isUpperCell
        isPositiveVertex
        isNegativeVertex
        isVertexInCell    
        areVertexCoherent        
        pathVertexes     
        isCellRight
        isCellLeft        
    end
    
    properties (Access = private)
        mesh
        orientation
        singularityCoord
    end
    
    methods (Access = public)
        
        function obj = CorrectorComputer(cParams)
            obj.init(cParams);
            obj.isUpperCell     = false(obj.mesh.nelem,1);
            obj.correctorValues = zeros(obj.mesh.nelem,obj.mesh.nnodeElem);
        end
                
        function phiV = compute(obj) 
            obj.computeCoherentOrientation();
            obj.computePathToBoundary();
            obj.createLeftRightPathElements();
            obj.computeReferenceCells();                        
            for ivertex = 1:length(obj.pathVertexes)
                obj.computeIsVertexInCell(ivertex);
                obj.isCellAnUpperCell(ivertex);                
                obj.computePositiveNegativeVertexes();
                obj.computeCorrector();
            end            
            phiV = obj.correctorValues;
        end
        
        function plot(obj)
            phi = obj.correctorValues;
            figure()
            s.mesh  = obj.mesh.createDiscontinousMesh();
            s.field = transpose(phi);
            n = NodalFieldPlotter(s);
            n.plot();
            shading interp
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.orientation        = cParams.orientation;
            obj.singularityCoord   = cParams.singularityCoord;
        end
        
        function computeCoherentOrientation(obj)
            s.mesh        = obj.mesh.createDiscontinuousMesh();
            s.orientation = obj.createDiscontinousField(obj.orientation);
            c = CoherentOrientationSelector(s);
            aC = c.isOrientationCoherent();
            obj.areVertexCoherent = aC;
        end                
        
        function computePathToBoundary(obj)
            s.mesh = obj.mesh;
            s.singularityCoord   = obj.singularityCoord;
            p = PathVertexesToBoundaryComputer(s);
            v = p.compute(); 
            obj.pathVertexes = v;
        end        

        function createLeftRightPathElements(obj)
            s.pathVertexes = obj.pathVertexes;
            s.mesh         = obj.mesh;
            l = LeftRightCellsOfPathToBoundaryComputer(s);
            [cR,cL] = l.compute();   
            l.plot();            
            obj.isCellLeft  = cL;
            obj.isCellRight = cR;
        end         
        
        function f = restrictToCell(obj,f)
            fV = obj.restrictToVertex(f);
            f  =  any(fV,2);
        end
        
        function fE = extendToAllVertexOfCell(obj,f)
            n = obj.mesh.nnodeElem;
            fE = repmat(f,1,n);
        end
        
        function computePositiveNegativeVertexes(obj)
            isCellPos  = obj.computePositiveCells();
            areVertPos = obj.extendToAllVertexOfCell(isCellPos);
            isP = obj.restrictToVertex(areVertPos);
            isN = obj.restrictToVertex(~areVertPos);
            obj.isPositiveVertex = isP;
            obj.isNegativeVertex = isN;
        end
        
        function isP = computePositiveCells(obj)
            isU = obj.isUpperCell;
            isRightDown = obj.isCellRight & ~isU;
            isLeftUpper = obj.isCellLeft  & isU;
            isP = isRightDown | isLeftUpper;
        end
        
        function computeCoherentAndNotCoherentVertex(obj)
            areC  = obj.areVertexCoherent;
            isCoh = obj.restrictToCell(areC);
            isNot = obj.restrictToCell(~areC);
            obj.isCoherent    = isCoh;
            obj.isNotCoherent = isNot;
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
            obj.hasCellSameCoherentOrientationAsReferenceCell(rCell);
            itHas    = obj.itHasSameCoherence;
            itHasNot = obj.itHasNotSameCoherence;
            isU = obj.isUpperCell;
            isU(itHas)    = isU(rCell);
            isU(itHasNot) = ~isU(rCell);
            obj.isUpperCell = isU;
        end
        
        function hasCellSameCoherentOrientationAsReferenceCell(obj,rCell)
            if obj.isCoherent(rCell)
                itHas    = obj.isCoherent;
                itHasNot = obj.isNotCoherent;
            else
                itHas    = obj.isNotCoherent;
                itHasNot = obj.isCoherent;
            end
            obj.itHasSameCoherence = itHas;
            obj.itHasNotSameCoherence = itHasNot;
        end
        
        function computeCorrector(obj)
            phi = obj.correctorValues;
            isP = obj.isPositiveVertex;
            isN = obj.isNegativeVertex;
            phi(isP) = 0.5;
            phi(isN) = -0.5;
            obj.correctorValues = phi;
        end
        
        function fV = restrictToVertex(obj,f)
            itIs = obj.isVertexInCell;
            fV = f & itIs;
        end        
        
        function fD = createDiscontinousField(obj,fValues)
%             s.connec = obj.mesh.connec;
%             s.type   = obj.mesh.type;
%             s.fNodes = fValues;
%             f = FeFunction(s);            
%             fD = f.computeDiscontinousField();
            s.fValues = fValues;
            s.connec = obj.mesh.connec;
            s.type   = obj.mesh.type;
            f = P1Function(s);
            s.mesh   = obj.mesh;
            s.connec = obj.mesh.connec;
            p = Projector_toP1Discontinuous(s);
            fD = p.project(f);
            fD = fD.fValues;
        end          
        
    end
    
 
    
end