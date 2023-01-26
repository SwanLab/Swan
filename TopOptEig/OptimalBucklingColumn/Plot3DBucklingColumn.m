classdef Plot3DBucklingColumn < handle
    
    properties (Access = private)
        backgroundMesh
        boundaryMesh
        levelSet
        unfittedMesh
        Rvec
    end
    
    properties (Access = private)
        radius
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
            obj.radius = cParams.radius;
            obj.coordinates    = cParams.coord;
        end
        
        function createBackgroundMesh(obj)
            nElem = length(obj.radius);
            x = linspace(0,2,20);
            y = linspace(0,2,20);
            z = linspace(0,4,nElem);
            [X,Y,Z] = meshgrid(x,y,z);   
            coord  = [X(:) Y(:) Z(:)];
            nnodes = size(coord,1);
            obj.computeRadiusVector(nElem,nnodes,coord);
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
            s.type       = 'cylinderBuckling';
            s.Radius     = obj.Rvec;
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

        function computeRadiusVector(obj,nElem,nnodes,coord)
            R = zeros(nnodes,1);
            k = size(coord,1)/nElem;
            c1 = 1;
            c2 = k;
            for i = 1:nElem
                R(c1:c2) = obj.radius(i);
                c1 = c1+k;
                c2 = c2+k;
            end
            obj.Rvec = R;
                
        end
        
    end
    
end