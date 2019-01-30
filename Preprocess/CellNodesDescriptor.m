classdef CellNodesDescriptor < handle
    
    properties (SetAccess = private)
        nodesInXmin
        nodesInXmax
        nodesInYmin
        nodesInYmax
        cornerNodes
    end
    
    properties (Access = private)
        x
        y
        xmin
        xmax
        ymin
        ymax
        allNodes
    end
    
    methods (Access = public)
        
        function obj = CellNodesDescriptor(coord)
            obj.init(coord);
            obj.obtainCellLimits();
            obj.obtainNodesInCellFaces();
            obj.obtainCornerNodes();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,coord)
            obj.x = coord(:,1);
            obj.y = coord(:,2);
            obj.allNodes(:,1) = 1:size(obj.x,1);
        end
        
        function obtainCellLimits(obj)
            obj.xmin = min(obj.x);
            obj.xmax = max(obj.x);
            obj.ymin = min(obj.y);
            obj.ymax = max(obj.y);
        end
        
        function obtainNodesInCellFaces(obj)
            obj.nodesInXmin = obj.getNodesInCellFace(obj.x,obj.xmin);
            obj.nodesInXmax = obj.getNodesInCellFace(obj.x,obj.xmax);
            obj.nodesInYmin = obj.getNodesInCellFace(obj.y,obj.ymin);
            obj.nodesInYmax = obj.getNodesInCellFace(obj.y,obj.ymax);
        end
        
        function obtainCornerNodes(obj)
            XminYmin = intersect(obj.nodesInXmin,obj.nodesInYmin);
            XminYmax = intersect(obj.nodesInXmin,obj.nodesInYmax);
            XmaxYmin = intersect(obj.nodesInXmax,obj.nodesInYmin);
            XmaxYmax = intersect(obj.nodesInXmax,obj.nodesInYmax);
            obj.cornerNodes = [XminYmin XminYmax XmaxYmin XmaxYmax]';
        end
        
        function nodesInFace = getNodesInCellFace(obj,x,xLim)
            dist2face = abs(x - xLim);
            isInFace = dist2face < 1e-13;
            nodesInFace = obj.allNodes(isInFace);
        end
    end
    
end

