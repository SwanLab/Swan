classdef STLExporter < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private) % Inputs
        unfittedMesh
    end
    
    properties (Access = private) % Calculated
        bgMesh, bdMesh
        jointMesh
        triMesh
        levelSet
    end
    
    methods (Access = public)
        
        function obj = STLExporter(cParams)
            obj.init(cParams)
        end

        function export(obj,filename)
            obj.joinInnerMeshes();
            obj.triangulateMesh();
            stlwrite(obj.triMesh, [filename '.stl'])
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.createSampleMeshes();
            obj.createSampleLevelSet();
            obj.createUnfittedMesh();
        end

        function joinInnerMeshes(obj)
            uM = obj.unfittedMesh;
            coordInner     = uM.innerMesh.mesh.coord;
            connecInner    = uM.innerMesh.mesh.connec;
            coordCutInner  = uM.innerCutMesh.mesh.coord;
            connecCutInner = uM.innerCutMesh.mesh.connec;
            ncoord = size(coordInner,1);
            connecCutInner = connecCutInner + ncoord;
            s.coord  = [coordInner;  coordCutInner];
            s.connec = [connecInner; connecCutInner];
            obj.jointMesh = Mesh(s);
        end

        function triangulateMesh(obj)
            P = obj.jointMesh.coord;
            T = obj.jointMesh.connec;
            obj.triMesh = triangulation(T,P);
        end

        function createSampleMeshes(obj)
            obj.createBackgroundMesh();
            obj.createBoundaryMesh();
        end

        function createBackgroundMesh(obj)
            x1 = linspace(-1,1,100);
            x2 = linspace(0,1,100);
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
        
    end
    
end