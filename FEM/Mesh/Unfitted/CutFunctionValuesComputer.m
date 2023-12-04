classdef CutFunctionValuesComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
      fValues
      newfValues
    end
    
    properties (Access = private)
        nodesInEdges
        isEdgeCut
    end
    
    properties (Access = private)
        backgroundfValues
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
            obj.backgroundfValues = cParams.oldfValues;
        end
        
        function computeCutPoints(obj)
            shapes = obj.computeShapes();
            shapeA = shapes(:,1);
            shapeB = shapes(:,2);
            node1 = obj.nodesInEdges(obj.isEdgeCut,1);
            node2 = obj.nodesInEdges(obj.isEdgeCut,2);
            fA  = obj.backgroundfValues(node1,:);
            fB  = obj.backgroundfValues(node2,:);
            fCut = zeros(size(fA));
            nnode = size(fA,2);
            for idim = 1:nnode
                fCut(:,idim) = fA(:,idim).*shapeA + fB(:,idim).*shapeB;
            end
            obj.newfValues = fCut;
        end
        
        function shapes = computeShapes(obj)
            m.type = 'LINE';
            m.coord  = [];
            m.connec = [];
            int = Interpolation.create(m,'LINEAR');
            xCutIso = obj.xCutEdgePoint';
            int.computeShapeDeriv(xCutIso);
            shapes = int.shape';
        end
        
        function computeCutMeshCoordinates(obj)
            obj.fValues = [obj.backgroundfValues;obj.newfValues];
        end
        
    end
    
end