classdef FreeFemFileReader < FileReader
    
    properties (Access = private)
        linesRead
    end
    
    methods (Access = public)
        
        function obj = FreeFemFileReader(path)
            obj.init(path)
        end
        
        function read(obj)
            obj.openFile();
            obj.readLines();
            obj.closeFile();
        end
        
        function d = getDataBase(obj)
            d = obj.linesRead;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,path)
            obj.filePath = path;
        end
        
        function readLines(obj)
            tline = fgetl(obj.fid);
            iline = 0;
            while ischar(tline)
                iline = iline + 1;
                obj.linesRead{iline} = tline;
                tline = fgetl(obj.fid);
            end
        end
    
    end

end