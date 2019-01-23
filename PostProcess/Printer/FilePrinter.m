classdef FilePrinter < handle
    
    properties (Access = protected, Abstract)
       ext 
    end
    
    properties (Access = protected)
        fileID
        fileName
        outFileName
        resultsDir        
        openMode = 'w'
        iter
        
        ndim
        etype
    end
    
    
    methods (Access = protected)
        
        function openFile(obj)
            obj.fileID = fopen(obj.fileName,obj.openMode);
        end
        
        function closeFile(obj)
            fclose(obj.fileID);
        end
        
        function createFileName(obj,iS) 
            fileNameWithNoPath = strcat(obj.outFileName,num2str(iS),'.flavia.',obj.ext);            
            obj.fileName = fullfile(obj.resultsDir,fileNameWithNoPath);
        end        
        
    end

    
end