classdef UnstructuredMeshCreator < handle
        
    properties (Access = private)
        xmin
        xmax
        ymin
        ymax
        nInteriorPoints
        nBoundaryPoints        
    end
    
    methods (Access = public)
        
        function obj = UnstructuredMeshCreator(cParams)
            obj.init(cParams)            
        end
        
        function m = create(obj)
            [xB,yB] = obj.createBoundaryVertex();
            [xI,yI] = obj.createInteriorVertex();            
            xT = [xB,xI];
            yT = [yB,yI];            
            s.coord(:,1) = xT;
            s.coord(:,2) = yT;
            s.connec = delaunay(s.coord);
            m = Mesh.create(s);
        end         
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.xmin = cParams.xmin;
            obj.xmax = cParams.xmax;
            obj.ymin = cParams.ymin;
            obj.ymax = cParams.ymax;
            obj.nInteriorPoints = cParams.nInteriorPoints;
            obj.nBoundaryPoints = cParams.nBoundaryPoints;            
        end
        
        function [xB,yB] = createBoundaryVertex(obj)
            [xD,yD] = obj.createDownVertex();
            [xU,yU] = obj.createUpperVertex();
            [xL,yL] = obj.createLeftVertex();
            [xR,yR] = obj.createRightVertex();
            xB = [xD,xU,xL,xR];
            yB = [yD,yU,yL,yR];
        end
        
        function [xI,yI] = createInteriorVertex(obj)
            n = obj.nInteriorPoints;
            xI = obj.createInteriorVertexCoord(obj.xmin,obj.xmax,n);
            yI = obj.createInteriorVertexCoord(obj.ymin,obj.ymax,n);
        end          
        
       function x = createRandVertexInClosedUnitInterval(obj,n)
            x0 = obj.createRandVertexInOpenUnitInterval(n-2);
            x  = [0,x0,1];
        end
        
        function x = createInteriorVertexCoord(obj,xMin,xMax,n)
            xUnit = obj.createRandVertexInOpenUnitInterval(n);
            x     = obj.scaleVertexCoord(xMin,xMax,xUnit);
        end
        
        function x = createBoundaryVertexCoord(obj,xMin,xMax,n)
            xUnit = obj.createRandVertexInClosedUnitInterval(n);
            x     = obj.scaleVertexCoord(xMin,xMax,xUnit);
        end
        
        function [x,y] = createDownVertex(obj)
            n = obj.nBoundaryPoints;
            x = obj.createBoundaryVertexCoord(obj.xmin,obj.xmax,n);            
            y = obj.ymin*ones(1,n);           
        end
        
        function [x,y] = createUpperVertex(obj)
            n = obj.nBoundaryPoints;
            x = obj.createBoundaryVertexCoord(obj.xmin,obj.xmax,n);            
            y = obj.ymax*ones(1,n);                  
        end
        
        function [x,y] = createLeftVertex(obj)
            n = obj.nBoundaryPoints;
            x = obj.xmin*ones(1,n);
            y = obj.createBoundaryVertexCoord(obj.ymin,obj.ymax,n);                        
        end
        
        function [x,y] = createRightVertex(obj)
            n = obj.nBoundaryPoints;
            x = obj.xmax*ones(1,n);                  
            y = obj.createBoundaryVertexCoord(obj.ymin,obj.ymax,n);                        
        end  
        
    end
    
    methods (Access = private, Static)
        
        function xV = createRandVertexInOpenUnitInterval(n)
            xV = rand(1,n);
        end        
        
        function x = scaleVertexCoord(xMin,xMax,xUnit)
            x = xMin + (xMax - xMin)*(xUnit);
        end
        
    end        
    
end