classdef DesignVariableStorerFromRes < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        levelSet
    end
    
    properties (Access = private)
        pathFile
        testName
        finalIter         
    end
    
    methods (Access = public)
        
        function obj = DesignVariableStorerFromRes(cParams)
            obj.init(cParams)            
        end
        
        function store(obj)
             for iter = 1:obj.finalIter
               obj.obtainLevelSet(iter);
               obj.saveLevelSet(iter);
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.finalIter = cParams.finalIter;
            obj.pathFile  = cParams.pathFile;
            obj.testName  = cParams.testName;
        end
        
        function obtainLevelSet(obj,iter)
            s.fileName = [obj.testName,num2str(iter)];
            s.folderPath = obj.pathFile;            
            wM = WrapperMshResFiles(s);
            wM.compute();
            obj.levelSet = wM.dataRes.DesignVar1;             
        end
        
        function saveLevelSet(obj,iter)
            x = obj.levelSet;
            outputName = ['/DesignVariable',num2str(iter)];
            outputFile = fullfile(obj.pathFile,outputName);
            save(outputFile,'x');
        end
        
    end
    
end