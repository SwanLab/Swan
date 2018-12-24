classdef testTopOptDensityPrinting < ...
            testNotShowingError & testTopOptComputation
    
    properties (Access = protected)
        testName = 'test_cantilever';  
        filesHaveChanged
        fileOutputName
        iter
    end
    
    methods (Access = public)
        
        function obj = testTopOptDensityPrinting()
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
            obj.fileOutputName = ['testTopOptDensityPrinting'];
            obj.iter = 0;
        end
        
        function compareFiles(obj)
            out   = obj.fileOutputName;
            name1 = ['tests/PrintingTests/PrintedFiles/testTopOptDensityPrinting_',num2str(obj.iter),'.flavia.res'];
            name2 = ['Output/',out,'/',out,'_',num2str(obj.iter),'.flavia.res'];
            command = ['diff ', name1, ' ', name2];
            [obj.filesHaveChanged,~] = system(command);
        end
        
        function print(obj)
            x         = obj.topOpt.x;
            optimizer = obj.topOpt.settings.optimizer;
            outName   = obj.fileOutputName;
            mesh      = obj.topOpt.mesh;
            results.density = x;
            results.iter = obj.iter;
            results.case_file = outName;
            postprocess = Postprocess_TopOpt.Create(optimizer);
            postprocess.print(mesh,results);
        end
        
    end
    
end

