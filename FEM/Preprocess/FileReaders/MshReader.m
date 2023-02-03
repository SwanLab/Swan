classdef MshReader < FileReader
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        strcontent
        contentCell
        coord
        connec
        nElem
        nNodes
        nNodeElem
        numbers
    end
    
    methods (Access = public)
        
        function obj = MshReader(cParams)
            obj.init(cParams)
        end

        function m = read(obj)
            obj.openFile();
            obj.splitIntoLines();
            obj.getDimensions();
            obj.extractNumbers();
            obj.getCoordAndConnec();
            obj.closeFile();
            m = obj.createMesh();
        end

        function d = getDataBase(obj)
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.filePath = cParams.filePath;
        end

        function splitIntoLines(obj)
            obj.strcontent = fileread(obj.filePath);
            obj.contentCell = splitlines(obj.strcontent);
        end

        function getDimensions(obj)
            fLCell = split(obj.contentCell{1});
            fLCoor = 3;
            lLCoor = find(strcmp(obj.contentCell, 'End Coordinates'))-1;
            fLConnec = find(strcmp(obj.contentCell, 'Elements')) +1;
            lLConnec  = length(obj.contentCell) - 2;
            obj.nNodes    = lLCoor - fLCoor + 1;
            obj.nElem     = lLConnec - fLConnec + 1;
            obj.nNodeElem = str2double(fLCell{end});
        end

        function extractNumbers(obj)
            expr = '(-?0?\.?\d*)';
            match = regexp(obj.strcontent,expr, 'tokens', 'dotexceptnewline');
            cellfunres2 = cellfun(@cell2mat,match,'UniformOutput',false);
            numbersAll = str2double(cellfunres2);
            obj.numbers = numbersAll(3:end);
        end

        function getCoordAndConnec(obj)
            nN  = obj.nNodes;
            nE  = obj.nElem;
            nNE = obj.nNodeElem;
            coordRaw  = reshape(obj.numbers(1:nN*4), [4, nN])';
            connecRaw = reshape(obj.numbers(nN*4+1:end), [nNE+1, nE])';
            if isequal(coordRaw(:,4),zeros(length(coordRaw),1))
                coordRaw = coordRaw(:,1:3);
            end
            obj.coord  = coordRaw(:,2:end);
            obj.connec = connecRaw(:,2:end);
        end
        
        function m = createMesh(obj)
            s.coord  = obj.coord;
            s.connec = obj.connec;
            m = Mesh(s);
        end
    end
    
end