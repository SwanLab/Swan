classdef CutFunctionValuesComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
      allValues
      cutValues
    end
    
    properties (Access = private)
        nodesInEdges
        isEdgeCut
    end
    
    properties (Access = private)
        fValues
        xCutEdgePoint
    end
    
    methods (Access = public)
        
        function obj = CutFunctionValuesComputer(cParams)
            obj.init(cParams)
        end
        
        function compute(obj)
            obj.computeCutPoints();
            obj.computeCutMeshCoordinates();
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.isEdgeCut       = cParams.isEdgeCut;
            obj.nodesInEdges    = cParams.nodesInEdges;
            obj.xCutEdgePoint   = cParams.xCutEdgePoint;
            obj.fValues         = cParams.fValues;
        end
        
        function computeCutPoints(obj)
            shapes = obj.computeShapes();
            shapeA = shapes(:,1);
            shapeB = shapes(:,2);
            node1 = obj.nodesInEdges(obj.isEdgeCut,1);
            node2 = obj.nodesInEdges(obj.isEdgeCut,2);
            fA  = obj.fValues(node1,:);
            fB  = obj.fValues(node2,:);
            fCut = zeros(size(fA));
            nnode = size(fA,2);
            for idim = 1:nnode
                fCut(:,idim) = fA(:,idim).*shapeA + fB(:,idim).*shapeB;
            end
            obj.cutValues = fCut;
        end
        
        function shapes = computeShapes(obj)
            type = 'LINE';
            int = Interpolation.create(type,'LINEAR');
            xCutIso = obj.xCutEdgePoint';
            shapes = int.computeShapeFunctions(xCutIso)';
        end
        
        function computeCutMeshCoordinates(obj)
            obj.allValues = [obj.fValues;obj.cutValues];
        end
        
    end
    
end