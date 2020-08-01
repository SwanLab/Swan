classdef CutCoordinatesComputer < handle
    
    properties (GetAccess = public, SetAccess = private)
      coord
      xCutPoints      
    end
    
    properties (Access = private)
        nodesInEdges
        isEdgeCut
    end
    
    properties (Access = private)
        backgroundCoord
        xCutEdgePoint
    end
    
    methods (Access = public)
        
        function obj = CutCoordinatesComputer(cParams)
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
            obj.backgroundCoord = cParams.coord;
        end
        
        function computeCutPoints(obj)
            shapes = obj.computeShapes();            
            shapeA = shapes(:,1);
            shapeB = shapes(:,2);            
            node1 = obj.nodesInEdges(obj.isEdgeCut,1);
            node2 = obj.nodesInEdges(obj.isEdgeCut,2);            
            xA  = obj.backgroundCoord(node1,:);
            xB  = obj.backgroundCoord(node2,:);
            xCut = zeros(size(xA));
            nnode = size(xA,2);
            for idim = 1:nnode
                xCut(:,idim) = xA(:,idim).*shapeA + xB(:,idim).*shapeB;
            end
            obj.xCutPoints = xCut;
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
            obj.coord = [obj.backgroundCoord;obj.xCutPoints];
        end
        
    end
    
end