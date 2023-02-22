classdef Plot3DBucklingColumn < handle
    
    properties (Access = private)
        backgroundMesh
        boundaryMesh
        levelSet
        unfittedMesh
    end
    
    properties (Access = private)
        value
        coordinates
    end
    
    methods (Access = public)
        
        function obj = Plot3DBucklingColumn(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            obj.createBackgroundMesh();
            obj.createBoundaryMesh();
            obj.createLevelSet();
            obj.createUnfittedMesh();
            obj.plotUnfittedMesh();
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.value = cParams.designVariableValue;
            obj.coordinates    = cParams.coord;
        end
        
        function createBackgroundMesh(obj)
            nElem = length(obj.coordinates)-1;
            x = linspace(0,3,20);
            y = linspace(0,3,20);
            z = linspace(0,4,nElem);
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
            s.type = 'rectangularColumn'; %'cylinderBuckling'/'holedCircle'/'rectangularColumn'
            s.desVarValue     = obj.value;
            s.coord      = obj.backgroundMesh.coord;
            s.ndim       = obj.backgroundMesh.ndim;
            lsCreator    = LevelSetCreator.create(s);
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