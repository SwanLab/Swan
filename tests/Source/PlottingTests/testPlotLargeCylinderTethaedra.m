classdef testPlotLargeCylinderTethaedra < testNotShowingError...
                                         & test
                                     
    properties (Access = private)
        backgroundMesh
        boundaryMesh
        levelSet
        unfittedMesh
    end
    
    methods (Access = public)
        
        function obj = testPlotLargeCylinderTethaedra() 
            obj.createBackgroundMesh();
            obj.createBoundaryMesh();
            obj.createLevelSet();
            obj.createUnfittedMesh();
            obj.plotUnfittedMesh();
        end
        
    end
    
    methods (Access = protected)        
       
        function hasPassed = hasPassed(obj)
            hasPassed = true;
        end        
        
    end
    
    methods (Access = private)

        function createBackgroundMesh(obj)
            x = linspace(0,1,10);
            y = linspace(0,1,10);
            z = linspace(0,2,20);
            [X,Y,Z] = meshgrid(x,y,z);   
            coord  = [X(:) Y(:) Z(:)];
            d = delaunayTriangulation(coord);
            s.connec = d.ConnectivityList;
            s.coord  = coord;
            obj.backgroundMesh = Mesh(s);            
        end
        
        function createBoundaryMesh(obj)
            s.backgroundMesh = obj.backgroundMesh;
            s.dimension = 1:obj.backgroundMesh.ndim;
            bC = BoundaryMeshCreatorFromRectangularBox(s);
            bM = bC.create();   
            obj.boundaryMesh = bM;
        end        
        
        function createLevelSet(obj)
            s.type = 'cylinder';
            s.fracRadius = 1.1;
            s.coord      = obj.backgroundMesh.coord;
            s.ndim       = obj.backgroundMesh.ndim;
            lsCreator = LevelSetCreator.create(s);
            obj.levelSet = lsCreator.getValue();            
        end        
        
        function createUnfittedMesh(obj)
            s.backgroundMesh = obj.backgroundMesh;
            s.boundaryMesh   = obj.boundaryMesh;
            uM = UnfittedMesh(s);            
            uM.compute(obj.levelSet);     
            obj.unfittedMesh = uM;
        end
        
        function plotUnfittedMesh(obj) 
            figure();
            obj.unfittedMesh.plot();
            view([1 1 1])            
        end
        
    end
end

