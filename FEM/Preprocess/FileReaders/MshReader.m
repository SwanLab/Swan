classdef MshReader < FileReader
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        contentCell
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = MshReader(cParams)
            obj.init(cParams)
        end

        function m = read(obj)
            obj.openFile();
            obj.splitIntoLines();
            obj.readCoords(); % 3 : find(strcmp(obj.contentCell, 'End Coordinates')) -1
            obj.readConnec(); % find(strcmp(obj.contentCell, 'Elements')) +1 : end-1
            obj.closeFile();
        end

        function d = getDataBase(obj)
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.filePath = cParams.filePath;
        end

        function splitIntoLines(obj)
            content = fileread(obj.filePath);
            obj.contentCell = splitlines(content);
        end

        function readCoords(obj)
            firstLine = 3;
            lastLine  = find(strcmp(obj.contentCell, 'End Coordinates'))-1;
        end

        function readConnec(obj)
            firstLine = find(strcmp(obj.contentCell, 'Elements')) +1;
            lastLine  = length(obj.contentCell) - 2;
            % {obj.contentCell{firstLine:lastLine}}'
        end
        
    end
    
end