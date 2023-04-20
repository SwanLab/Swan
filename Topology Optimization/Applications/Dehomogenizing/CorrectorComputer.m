classdef CorrectorComputer < handle
       
    properties (Access = private)
        upperCellFunction
        correctorFunction
        correctorValues
        

        isUpperCell
        isPositiveVertex
        isNegativeVertex
        isVertexInCell    
        pathVertexes     
        isCellRight
        isCellLeft
        singularCoord
    end
    
    properties (Access = private)
        mesh
        singularElement
        areCoherent
    end
    
    methods (Access = public)
        
        function obj = CorrectorComputer(cParams)
            obj.init(cParams);
            obj.correctorValues = zeros(obj.mesh.nelem,obj.mesh.nnodeElem);
        end
                
        function cF = compute(obj) 
            obj.computeSingularCoord();
            obj.computePathToBoundary();
            obj.computeLeftRightPathElements();
            obj.computeUpperLowerCell();
            for ivertex = 1:length(obj.pathVertexes)
                obj.computeIsVertexInCell(ivertex);
                obj.computePositiveNegativeVertexes();
                obj.computeCorrector();
            end    
            obj.createCorrectorFunction();
            cF = obj.correctorFunction;            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh            = cParams.mesh;
            obj.areCoherent     = cParams.isCoherent;
            obj.singularElement = cParams.singularElement;
        end

        function computeSingularCoord(obj)
            isS = obj.singularElement;
            sC  = obj.mesh.computeBaricenter();
            sC  = transpose(sC);
            sC  = sC(isS,:); 
            obj.singularCoord = sC;
        end
        
         function computePathToBoundary(obj)
            s.mesh = obj.mesh;
            s.singularityCoord   = obj.singularCoord;
            p = PathVertexesToBoundaryComputer(s);
            v = p.compute(); 
            obj.pathVertexes = v(1:end);
        end        

        function computeLeftRightPathElements(obj)
            s.pathVertexes = obj.pathVertexes;
            s.mesh         = obj.mesh;
            s.singularElement = obj.singularElement;
            l = LeftRightCellsOfPathToBoundaryComputer(s);
            [cR,cL] = l.compute();   
          %  l.plot();            
            obj.isCellLeft  = cL;
            obj.isCellRight = cR;
        end    

        function computeUpperLowerCell(obj)
            s.pathVertexes    = obj.pathVertexes;
            s.mesh            = obj.mesh;
            s.areCoherent     = obj.areCoherent;            
            u = UpperLowerCellComputer(s);
            u.compute();
            obj.isUpperCell = u.isUpperCell;
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
            isU = squeeze(obj.isUpperCell.fValues);
            isRightDown = obj.isCellRight & ~isU;
            isLeftUpper = obj.isCellLeft  & isU;
            isP = isRightDown | isLeftUpper;
        end
         
        function computeIsVertexInCell(obj,ivertex)
            vI = obj.pathVertexes(ivertex);
            vertexInCell  = obj.mesh.connec;
            itIs = (vertexInCell == vI);
            obj.isVertexInCell = itIs;
        end
        
        function computeCorrector(obj)
            phi = obj.correctorValues;
            isP = obj.isPositiveVertex;
            isN = obj.isNegativeVertex;
            phi(isP) = 0.5;
            phi(isN) = -0.5;
            obj.correctorValues = phi;
        end

        function createCorrectorFunction(obj)
            s.fValues = permute(obj.correctorValues, [3, 2, 1]);
            s.mesh    = obj.mesh;
            f = P1DiscontinuousFunction(s);
            obj.correctorFunction = f;
        end    
        
        function fV = restrictToVertex(obj,f)
            itIs = obj.isVertexInCell;
            fV = f & itIs;
        end                
        
    end
    
end