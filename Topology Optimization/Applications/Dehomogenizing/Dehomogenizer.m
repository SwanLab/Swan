classdef Dehomogenizer < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        boundaryMesh
        uMesh        
    end
    
    properties (Access = private)
        cellLevelSetParams
        backgroundMesh
        nCells
        theta
        mesh
    end
    
    methods (Access = public)
        
        function obj = Dehomogenizer(cParams)
            obj.init(cParams)            
        end
        
        function compute(obj)
            obj.createBoundaryMesh();                       
            obj.createUnfittedMesh();            
        end
        
        function plot(obj)
            obj.plotOrientation();
            obj.plotStructure();
            obj.plotComponents();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.backgroundMesh     = cParams.backgroundMesh;
            obj.nCells             = cParams.nCells;
            obj.theta              = cParams.theta;
            obj.cellLevelSetParams = cParams.cellLevelSetParams;            
            obj.mesh               = cParams.mesh;
        end
        
        function createBoundaryMesh(obj)
            sB.backgroundMesh = obj.backgroundMesh;
            sB.dimension = 1:3;
            sB.type = 'FromReactangularBox';
            bMc = BoundaryMeshCreator.create(sB);
            obj.boundaryMesh  = bMc.create();            
        end
        
        function createUnfittedMesh(obj)
            ls = obj.createLevelSet();                                    
            s.boundaryMesh   = obj.boundaryMesh;
            s.backgroundMesh = obj.backgroundMesh;
            obj.uMesh = UnfittedMesh(s);
            obj.uMesh.compute(ls);
        end               
        
        function plotOrientation(obj)
            figure()
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            t  = obj.theta;
            ct = cos(t(:,1));
            st = sin(t(:,1));
            quiver(x,y,ct,st)
        end
        
        function plotStructure(obj)
            figure()
            obj.uMesh.plotStructureInColor('black');
        end
        
        function plotComponents(obj)
            s.unfittedMesh = obj.uMesh;
            sp = UnfittedMeshSplitter(s);
            sp.split();                        
        end

        function ls = createLevelSet(obj)
            s.coord = obj.backgroundMesh.coord;            
            s.type   = 'periodicAndOriented';            
            s.backgroundMesh   = obj.backgroundMesh;
            s.mesh   = obj.mesh;
            s.ndim   = 2;            
            s.angle  = obj.theta;
            s.nCells = obj.nCells;
            s.cellLevelSetParams = obj.cellLevelSetParams;
            levelSet = LevelSetCreator.create(s);            
            ls = levelSet.getValue();   
        end             
        
    end
    
end