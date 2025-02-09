classdef FileReader < handle
    
    properties (Access = protected)
        filePath
        fid
    end
    
    
    methods (Access = protected)
        
        function openFile(obj)
            obj.fid = fopen(obj.filePath);
        end
        
        function closeFile(obj)
            fclose(obj.fid);
        end
        
    end
    
    methods (Abstract, Access = public)
        read(obj)
        getDataBase(obj)
    end
    
end