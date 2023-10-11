classdef ConnectivityComputer < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
       mesh
       levelSet
       uMesh
       boundaryMesh
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = ConnectivityComputer()
            obj.init();
            obj.createMesh();
            obj.createBoundaryMesh();
            obj.createLevelSet();
            obj.createUnfittedMesh();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            
        end
        
        function createMesh(obj)
            x1 = linspace(-1,1,20);
            x2 = linspace(0,1,20);
            [xv,yv] = meshgrid(x1,x2);
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            m = Mesh(s);
            figure()
            m.plot();
            obj.mesh = m;
        end

        function createBoundaryMesh(obj)
            s.backgroundMesh = obj.mesh;
            s.dimension      = 1:3;
            s.type = 'FromReactangularBox';
            bM = BoundaryMeshCreator.create(s);
            m  = bM.create();
            obj.boundaryMesh = m;
        end

        function createLevelSet(obj)
            s.type       = 'circleInclusion';
            s.mesh       = obj.mesh;
            s.ndim       = 2;
            s.fracRadius = 0.4;
            s.coord      = obj.mesh.coord;
            ls = LevelSetCreator.create(s);
            obj.levelSet = ls.getValue();
        end
        
        function createUnfittedMesh(obj)
            s.boundaryMesh   = obj.boundaryMesh;
            s.backgroundMesh = obj.mesh;
            m = UnfittedMesh(s);
            m.compute(obj.levelSet);
            figure
            m.plot()
            obj.uMesh = m;
        end
        
    end
    
end