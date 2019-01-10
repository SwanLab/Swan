classdef testPrintingDescriptor < handle
    
    properties (Access = protected, Abstract)
        fileOutputName
        iter
    end
    
    properties (Access = protected)
       filesHaveChanged 
    end
    
    methods (Access = protected)
        
        function compareFiles(obj)
            hasMshChanged = obj.compareFile('.msh');
            hasResChanged = obj.compareFile('.res');
            obj.filesHaveChanged = hasMshChanged || hasResChanged;
        end
        
        function hasChanged = compareFile(obj,extension)
            out   = obj.fileOutputName;
            fullOutSavedName = [out,'.flavia',extension];
            fullOutName = [out,num2str(obj.iter),'.flavia',extension];
            savedPrintedFile = fullfile('tests','PrintingTests','PrintedFiles',fullOutSavedName);
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