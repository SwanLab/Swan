classdef DesignVariableAndMeshSaverFromRes < handle
    
    properties (Access = private)
        levelSet
        backgroundMesh
        iterationsNumber
    end
    
    properties (Access = private)
        testPath
        testName
    end
    
    methods (Access = public)
        
        function obj = DesignVariableAndMeshSaverFromRes(cParams)
            obj.init(cParams)
        end
        
        function save(obj)
            obj.readIterNumbers()
            for iter = 1:length(obj.iterationsNumber)                
                obj.obtainLevelSetAndMesh(iter);
                obj.saveLevelSetAndMesh(iter);
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.testPath  = cParams.testPath;
            obj.testName  = cParams.testName;
        end
        
        function readIterNumbers(obj)
            files = dir(fullfile(obj.testPath,'*.res'));
            nFiles = numel(files);
            iterN  = zeros(nFiles,1);
            for ifiles = 1:nFiles
                file = files(ifiles);
                iter = obj.getIter(file);
                iterN(ifiles) = iter;
            end
            obj.iterationsNumber = sort(iterN);
        end
        
        
        function obtainLevelSetAndMesh(obj,iter)
            it = obj.iterationsNumber(iter);            
            s.fileName = [obj.testName,num2str(it)];
            s.folderPath = fullfile(obj.testPath);
            wM = WrapperMshResFiles(s);
            wM.compute();
            obj.levelSet       = wM.dataRes.DesignVar1;
            obj.backgroundMesh = wM.mesh;
        end
        
        function saveLevelSetAndMesh(obj,iter)
            x = obj.levelSet;
            mesh = obj.backgroundMesh;
            outputName = ['DesignVariable',num2str(iter)];
            outputFile = fullfile(obj.testPath,outputName);
            save(outputFile,'x','mesh');
        end
        
    end
    
    methods (Access = private, Static)
        
        function iter = getIter(file)
            fName = file.name;
            iterCell = regexp(fName,'[0-9]','match');
            iterM    = cell2mat(iterCell);
            iter     = str2double(iterM);
        end
        
    end
    
end