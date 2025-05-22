classdef MshReader < FileReader
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        contentCell
        coord
        connec
        nElem
        nNodes
        nNodeElem
        fLCoor, lLCoor
        fLConnec, lLConnec
    end
    
    methods (Access = public)
        
        function obj = MshReader(cParams)
            obj.init(cParams)
        end

        function m = read(obj)
            obj.openFile();
            obj.splitIntoLines();
            obj.closeFile();
            obj.getDimensions();
            obj.getCoordinates();
            obj.getConnectivity();
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
            strcontent = fileread(obj.filePath);
            obj.contentCell = splitlines(strcontent);
        end

        function getDimensions(obj)
            fLCell = split(obj.contentCell{1});
            obj.fLCoor   = 3;
            obj.lLCoor   = find(strcmp(obj.contentCell, 'End Coordinates'))-1;
            obj.fLConnec = find(strcmp(obj.contentCell, 'Elements')) +1;
            obj.lLConnec = length(obj.contentCell) - 2;
            obj.nNodes    = obj.lLCoor - obj.fLCoor + 1;
            obj.nElem     = obj.lLConnec - obj.fLConnec + 1;
            obj.nNodeElem = str2double(fLCell{end});
        end
        
        function getCoordinates(obj)
            strCoord = sprintf('%s ', obj.contentCell{obj.fLCoor:obj.lLCoor});
            fmtCoord = '%d %f %f %f';
            scanCoord = sscanf(strCoord, fmtCoord);
            coords = reshape(scanCoord, [4, obj.nNodes])';
            if isequal(coords(:,4),zeros(length(coords),1))
                coords = coords(:,1:3);
            end
            obj.coord = coords(:,2:end);
        end

        function getConnectivity(obj)
            fmtConnec = repmat('%d ', [1,obj.nNodeElem]);
            strConnec = sprintf('%s ', obj.contentCell{obj.fLConnec:obj.lLConnec});
            scanConnec = sscanf(strConnec, fmtConnec);
            conn = reshape(scanConnec, [obj.nNodeElem+1, obj.nElem-1])';
            obj.connec = conn(:,2:end);
        end

        function m = createMesh(obj)
            s.coord  = obj.coord;
            s.connec = obj.connec;
            m = Mesh.create(s);
        end
    end
    
end