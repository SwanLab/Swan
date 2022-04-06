classdef FilePrinter < handle
    
    properties (Access = protected)
        fileID
        fileName
        openMode = 'w'
    end
    
    methods (Access = protected)
        
        function openFile(obj)
            obj.fileID = fopen(obj.fileName,obj.openMode);
        end
        
        function closeFile(obj)
            fclose(obj.fileID);
        end
        
    end
    
end