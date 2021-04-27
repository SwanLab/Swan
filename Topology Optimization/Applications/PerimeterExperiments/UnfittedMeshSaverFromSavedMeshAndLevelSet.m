classdef UnfittedMeshSaverFromSavedMeshAndLevelSet < handle
    
    properties (Access = private)
        levelSet
        backgroundMesh
        boundaryMesh
        unfittedMesh
    end
    
    properties (Access = private)
        testPath
    end
    
    methods (Access = public)
        
        function obj = UnfittedMeshSaverFromSavedMeshAndLevelSet(cParams)
            obj.init(cParams)
            
        end
        
        function save(obj)
            nIter = obj.computeNiter();
            for iter = 1:nIter
                obj.loadDesignVariableAndMesh(iter);
                obj.createBoundaryMesh();
                obj.createUnfittedMesh(); 
                obj.saveUnfittedMesh(iter);
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.testPath  = cParams.testPath;
        end
        
        function nFiles = computeNiter(obj)
            files = dir(fullfile(obj.testPath,'DesignVariable*.mat'));
            nFiles = numel(files);            
        end
        
        function loadDesignVariableAndMesh(obj,iter)
            fullPath = fullfile(obj.testPath,'DesignVariable');
            d = load([fullPath,num2str(iter),'.mat']);
            obj.levelSet        = d.x;
            obj.backgroundMesh  = d.mesh;
        end
        
        function createBoundaryMesh(obj)
            sB.backgroundMesh = obj.backgroundMesh;
            sB.dimension = 1:3;
            sB.type = 'FromReactangularBox';
            bMc = BoundaryMeshCreator.create(sB);
            obj.boundaryMesh  = bMc.create();
        end                
        
        function createUnfittedMesh(obj)
            s.backgroundMesh  = obj.backgroundMesh;
            s.boundaryMesh    = obj.boundaryMesh;
            u = UnfittedMesh(s);
            u.compute(obj.levelSet);
            obj.unfittedMesh = u;
        end

        function saveUnfittedMesh(obj,iter)
            u = obj.unfittedMesh;
            outputName = ['UnfittedMesh',num2str(iter)];
            outputFile = fullfile(obj.testPath,outputName);            
            save(outputFile,'u');
        end
        
    end
    
end