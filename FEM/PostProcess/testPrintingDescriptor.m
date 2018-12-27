classdef testPrintingDescriptor < handle
    
    properties (Access = protected, Abstract)
        filesHaveChanged 
        fileOutputName
        iter
    end
    
    methods (Access = protected)
        
        function compareFiles(obj)
            hasMshChanged = obj.compareFile('.msh');
            hasResChanged = obj.compareFile('.res');
            obj.filesHaveChanged = hasMshChanged || hasResChanged;
        end
        
        function hasChanged = compareFile(obj,extension)
            out   = obj.fileOutputName;
            fullOutName = [out,num2str(obj.iter),'.flavia',extension];
            savedPrintedFile = fullfile('tests','PrintingTests','PrintedFiles',fullOutName);
            outputFile = fullfile('Output',out,fullOutName);
            command = ['diff ', savedPrintedFile, ' ', outputFile];
            [hasChanged,~] = system(command);
        end
        
        function hasPassed = hasPassed(obj)
            hasPassed = ~obj.filesHaveChanged;
        end
        
        function selectComputedVar(obj)
        end
        
    end

    
end    