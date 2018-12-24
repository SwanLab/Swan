classdef testTopOptLevelSetPrinting < ...
            testNotShowingError & testTopOptComputation
    
    properties (Access = protected)
        testName = 'test_gripping';  
        filesHaveChanged
        fileOutputName
        iter
    end
    
    methods (Access = public)
        
        function obj = testTopOptLevelSetPrinting()
            obj.init()
            obj.print()
            obj.compareFiles()
        end
    end
    
    methods (Access = protected)
        
        function hasPassed = hasPassed(obj)
            hasPassed = ~obj.filesHaveChanged;
        end
        
        function selectComputedVar(obj)
        end
        
    end

    methods (Access = private)
        
        function init(obj)
            obj.fileOutputName = ['testTopOptLevelSetPrinting'];
            obj.iter = 0;
        end
        
        function compareFiles(obj)
            hasMshChanged = obj.compareFile('.msh');
            hasResChanged = obj.compareFile('.res');
            obj.filesHaveChanged = hasMshChanged || hasResChanged;
        end
        
        function hasChanged = compareFile(obj,extension)
            out   = obj.fileOutputName;
            name1 = ['tests/PrintingTests/PrintedFiles/testTopOptLevelSetPrinting_',num2str(obj.iter),'.flavia',extension];
            name2 = ['Output/',out,'/',out,'_',num2str(obj.iter),'.flavia',extension];
            command = ['diff ', name1, ' ', name2];
            [hasChanged,~] = system(command);
        end
        
        function print(obj)
            x         = obj.topOpt.x;
            optimizer = obj.topOpt.settings.optimizer;
            outName   = obj.fileOutputName;
            mesh      = obj.topOpt.mesh;
            results.design_variable = x;
            results.iter = obj.iter;
            results.case_file = outName;
            postprocess = Postprocess_TopOpt.Create(optimizer);
            postprocess.print(mesh,results);
        end
        
    end
    
end

