classdef test_InnerMeshExporter < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        bgMesh, bdMesh
        levelSet
        unfittedMesh
    end
    
    methods (Access = public)
        
        function obj = test_InnerMeshExporter()
            obj.createSampleMeshes();
            obj.createSampleLevelSet();
            obj.createUnfittedMesh();
%             obj.printUnfittedMesh();
%             obj.exportMshGiD()
            obj.exportUsingExporter();
        end
        
    end
    
    methods (Access = private)
        
        function createSampleMeshes(obj)
            obj.createBackgroundMesh();
            obj.createBoundaryMesh();
        end

        function createBackgroundMesh(obj)
            x1 = linspace(-1,1,70);
            x2 = linspace(0,1,70);
            [xv,yv] = meshgrid(x1,x2);
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            obj.bgMesh = Mesh(s);
        end

        function createBoundaryMesh(obj)
            sB.backgroundMesh = obj.bgMesh;
            sB.dimension      = 1:3;
            sB.type = 'FromReactangularBox';
            bMc = BoundaryMeshCreator.create(sB);
            obj.bdMesh  = bMc.create();
        end

        function createSampleLevelSet(obj)
            s.type       = 'circleInclusion';
            s.mesh       = obj.bgMesh;
            s.ndim       = 2;
            s.fracRadius = 0.4;
            s.coord      = obj.bgMesh.coord;
            ls = LevelSetCreator.create(s);
            obj.levelSet = ls.getValue();
        end

        function createUnfittedMesh(obj)
            s.boundaryMesh   = obj.bdMesh;
            s.backgroundMesh = obj.bgMesh;
            uMesh = UnfittedMesh(s);
            uMesh.compute(obj.levelSet);
            obj.unfittedMesh = uMesh;
%             uMesh.plot();
        end

        function printUnfittedMesh(obj)
            filename = 'sample_test_innermeshexporter';
            obj.unfittedMesh.print(filename);
        end

        function exportMshGiD(obj)
            s.filename        = 'hellothere';
            s.gidProjectPath  = '/home/ton/test_micro_project.gid';
            s.meshElementSize = '0.0707107';
            s.meshFileName    = 'hmmmm22';
            s.swanPath        = '/home/ton/Github/Swan/';
            s.gidPath         = '/home/ton/GiDx64/gid-16.1.2d/';
            obj.unfittedMesh.exportGiD(s);
        end
        
        function exportUsingExporter(obj)
            m = obj.unfittedMesh.exportInnerMesh();
            m.plot();
        end

    end
    
end