classdef TestNaca < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        refMesh
        levelSet
        mesh
        meshGoodCond
    end
    
    methods (Access = public)
        
        function obj = TestNaca()
            obj.createReferenceMesh();
            obj.createLevelSet();
            obj.createFluidMesh();
            %obj.createFluidMeshGoodConditioning(); Maybe in the future if necessary
        end
        
    end
    
    methods (Access = private)
                
        function createReferenceMesh(obj)      
            obj.refMesh = TriangleMesh(2,1,150,75);
        end
        
        function createLevelSet(obj)
            s.type = 'Naca';
            s.xLE  = 0.5;
            s.yLE  = 0.5;

            s.chord = 1;
            s.p     = 0.5;
            s.m     = 0.02;
            s.t     = 0.12;


            g            = GeometricalFunction(s);
            lsFun        = g.computeLevelSetFunction(obj.refMesh);
            obj.levelSet = lsFun.fValues;
        end
        
        function createFluidMesh(obj)
            s.backgroundMesh = obj.refMesh;
            s.boundaryMesh   = obj.refMesh.createBoundaryMesh();
            uMesh            = UnfittedMesh(s);
            uMesh.compute(obj.levelSet);       
            obj.mesh = uMesh.createInnerMesh();
        end

        function createFluidMeshGoodConditioning(obj)
            points = obj.mesh.coord;
            T      = alphaShape(points,.007);
            DT     = alphaTriangulation(T);

            s.connec         = DT;
            s.coord          = points;
            obj.meshGoodCond = Mesh.create(s);
        end

    end

end