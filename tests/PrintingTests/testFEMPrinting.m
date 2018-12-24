classdef testFEMPrinting < ...
            testNotShowingError & testFemComputation
    
    properties (Access = protected)
        testName = 'test2d_quad';  
        filesHaveChanged
        fileOutputName
    end
    
    methods (Access = public)
        
        function obj = testFEMPrinting()
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
            obj.fileOutputName = ['testFemPrinting'];
        end
        
        function compareFiles(obj)
            hasMshChanged = obj.compareFile('.msh');
            hasResChanged = obj.compareFile('.res');
            obj.filesHaveChanged = hasMshChanged || hasResChanged;
        end
        
        function hasChanged = compareFile(obj,extension)
            out   = obj.fileOutputName;
            name1 = ['tests/PrintingTests/PrintedFiles/testFemPrinting_u_1.flavia',extension];
            name2 = ['Output/',out,'/',out,'_u_1.flavia',extension];
            command = ['diff ', name1, ' ', name2];
            [hasChanged,~] = system(command);
        end        
        
        
        function print(obj)
            var     = obj.fem.variables;
            outName = obj.fileOutputName;
            postprocess = Postprocess_PhysicalProblem;
            postprocess.print(obj.fem,outName,var);
        end
        
    end
    
end

