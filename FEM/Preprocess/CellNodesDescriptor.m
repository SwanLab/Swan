classdef CellNodesDescriptor < handle
    
    properties (SetAccess = private)
        nodesInXmin
        nodesInXmax
        nodesInYmin
        nodesInYmax
        nodesInZmin
        nodesInZmax
        cornerNodes
    end
    
    properties (Access = private)
        x
        y
        z
        xmin
        xmax
        ymin
        ymax
        zmin
        zmax
        allNodes
    end
    
    methods (Access = public)
        
        function obj = CellNodesDescriptor(coord)
            switch size(coord,2)
                case 2
                    obj.init2D(coord);
                    obj.obtainCellLimits2D();
                    obj.obtainNodesInCellFaces2D();
                    obj.obtainCornerNodes2D();
                case 3
                    obj.init3D(coord);
                    obj.obtainCellLimits3D();
                    obj.obtainNodesInCellFaces3D();
                    obj.obtainCornerNodes3D();
            end
        end
        
    end
    
    methods (Access = private)
        
        function init2D(obj,coord)
            obj.x = coord(:,1);
            obj.y = coord(:,2);
            obj.allNodes(:,1) = 1:size(obj.x,1);
        end
        
        function init3D(obj,coord)
            obj.x = coord(:,1);
            obj.y = coord(:,2);
            obj.z = coord(:,3);
            obj.allNodes(:,1) = 1:size(obj.x,1);
        end
        
        function obtainCellLimits2D(obj)
            obj.xmin = min(obj.x);
            obj.xmax = max(obj.x);
            obj.ymin = min(obj.y);
            obj.ymax = max(obj.y);
        end
        
        function obtainCellLimits3D(obj)
            obj.xmin = min(obj.x);
            obj.xmax = max(obj.x);
            obj.ymin = min(obj.y);
            obj.ymax = max(obj.y);
            obj.zmin = min(obj.z);
            obj.zmax = max(obj.z);
        end
        
        function obtainNodesInCellFaces2D(obj)
            obj.nodesInXmin = obj.getNodesInCellFace(obj.x,obj.xmin);
            obj.nodesInXmax = obj.getNodesInCellFace(obj.x,obj.xmax);
            obj.nodesInYmin = obj.getNodesInCellFace(obj.y,obj.ymin);
            obj.nodesInYmax = obj.getNodesInCellFace(obj.y,obj.ymax);
        end
        
        function obtainNodesInCellFaces3D(obj)
            obj.nodesInXmin = obj.getNodesInCellFace(obj.x,obj.xmin);
            obj.nodesInXmax = obj.getNodesInCellFace(obj.x,obj.xmax);
            obj.nodesInYmin = obj.getNodesInCellFace(obj.y,obj.ymin);
            obj.nodesInYmax = obj.getNodesInCellFace(obj.y,obj.ymax);
            obj.nodesInZmin = obj.getNodesInCellFace(obj.z,obj.zmin);
            obj.nodesInZmax = obj.getNodesInCellFace(obj.z,obj.zmax);
        end
        
        function obtainCornerNodes2D(obj)
            XminYmin = intersect(obj.nodesInXmin,obj.nodesInYmin);
            XminYmax = intersect(obj.nodesInXmin,obj.nodesInYmax);
            XmaxYmin = intersect(obj.nodesInXmax,obj.nodesInYmin);
            XmaxYmax = intersect(obj.nodesInXmax,obj.nodesInYmax);
            obj.cornerNodes = [XminYmin XminYmax XmaxYmin XmaxYmax]';
        end
        
        function obtainCornerNodes3D(obj)
            XminYmin = intersect(obj.nodesInXmin,obj.nodesInYmin);
            XminYmax = intersect(obj.nodesInXmin,obj.nodesInYmax);
            XmaxYmin = intersect(obj.nodesInXmax,obj.nodesInYmin);
            XmaxYmax = intersect(obj.nodesInXmax,obj.nodesInYmax);

            XminYminZmin = intersect(XminYmin,obj.nodesInZmin);
            XminYminZmax = intersect(XminYmin,obj.nodesInZmax);

            XminYmaxZmin = intersect(XminYmax,obj.nodesInZmin);
            XminYmaxZmax = intersect(XminYmax,obj.nodesInZmax);

            XmaxYminZmin = intersect(XmaxYmin,obj.nodesInZmin);
            XmaxYminZmax = intersect(XmaxYmin,obj.nodesInZmax);

            XmaxYmaxZmin = intersect(XmaxYmax,obj.nodesInZmin);
            XmaxYmaxZmax = intersect(XmaxYmax,obj.nodesInZmax);

            obj.cornerNodes = [XminYminZmin XminYminZmax XminYmaxZmin ...
                XminYmaxZmax XmaxYminZmin XmaxYminZmax XmaxYmaxZmin ...
                XmaxYmaxZmax]';
        end
        
        function nodesInFace = getNodesInCellFace(obj,x,xLim)
            dist2face = abs(x - xLim);
            isInFace = dist2face < 1e-13;
            nodesInFace = obj.allNodes(isInFace);
        end
    end
    
end

