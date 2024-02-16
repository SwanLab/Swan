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
            s.xmin = obj.xmin;
            s.xmax = obj.xmax;
            s.ymin = obj.ymin;
            s.ymax = obj.ymax;
            s.nInteriorPoints = obj.nInteriorPoints;
            s.nBoundaryPoints = obj.nBoundaryPoints;
            uM = UnstructuredMeshCreator(s);
            m = uM.create();
            obj.mesh = m;
        end
        
        function createBenchMarkMesh(obj)
            xv = linspace(obj.xmin,obj.xmax,3);
            yv = linspace(obj.ymin,obj.ymax,4);
            [X,Y] = meshgrid(xv,yv);
            [F,V] = mesh2tri(X,Y,zeros(size(X)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            m = Mesh.create(s);
            obj.mesh = m;
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

    
end