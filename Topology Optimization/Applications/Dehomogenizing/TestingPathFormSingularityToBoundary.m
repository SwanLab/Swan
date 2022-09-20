classdef TestingPathFormSingularityToBoundary < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        xmin
        xmax
        ymin
        ymax
        nInteriorPoints
        nBoundaryPoints
    end
    
    properties (Access = private)
        mesh
        singularityCoord
        boundaryPointCoord
        pathVertexes
        cellRight
        cellLeft        
    end
    
    methods (Access = public)
        
        function obj = TestingPathFormSingularityToBoundary()
            obj.init();
            obj.createMesh();
            %obj.createBenchMarkMesh()
            obj.createSingularityPoint();
            obj.createFinalPathPoint();
            obj.createPathToBoundary();
            obj.createLeftRightPathElements();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.xmin = 0;
            obj.xmax = 2;
            obj.ymin = 0;
            obj.ymax = 3;
            obj.nInteriorPoints = 20000;
            obj.nBoundaryPoints = 100;
        end
        
        function createMesh(obj)
            [xB,yB] = obj.createBoundaryVertex();
            [xI,yI] = obj.createInteriorVertex();            
            xT = [xB,xI];
            yT = [yB,yI];            
            s.coord(:,1) = xT;
            s.coord(:,2) = yT;
            s.connec = delaunay(s.coord);
            m = Mesh(s);
            obj.mesh = m;
        end        
        
        function createBenchMarkMesh(obj)
            xv = linspace(obj.xmin,obj.xmax,3);
            yv = linspace(obj.ymin,obj.ymax,4);
            [X,Y] = meshgrid(xv,yv);
            [F,V] = mesh2tri(X,Y,zeros(size(X)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            m = Mesh(s);
            obj.mesh = m;
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
        
        function [xI,yI] = createInteriorVertex(obj)
            n = obj.nInteriorPoints;
            xI = obj.createInteriorVertexCoord(obj.xmin,obj.xmax,n);
            yI = obj.createInteriorVertexCoord(obj.ymin,obj.ymax,n);
        end
        
        function [xB,yB] = createBoundaryVertex(obj)
            [xD,yD] = obj.createDownVertex();
            [xU,yU] = obj.createUpperVertex();
            [xL,yL] = obj.createLeftVertex();
            [xR,yR] = obj.createRightVertex();
            xB = [xD,xU,xL,xR];
            yB = [yD,yU,yL,yR];
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
        
        function createSingularityPoint(obj)
            obj.singularityCoord = [0.59 2.7];
        end
        
        function createFinalPathPoint(obj)
            obj.boundaryPointCoord = [2 1];
            obj.boundaryPointCoord = [0.59 0];            
        end        
        
        function createPathToBoundary(obj)
            s.mesh = obj.mesh;
            s.singularityCoord   = obj.singularityCoord;
            s.boundaryPointCoord = obj.boundaryPointCoord;
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
            obj.cellLeft  = cL;
            obj.cellRight = cR;            
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