classdef WrapperMshFile < FileReader

    properties (Access = private)
        dataBase
        lineNumber0
        tScan
        dimension
        nnode
    end

    properties (Access = private)
        lineNumber
    end

    methods (Access = public)

        function obj = WrapperMshFile(cParams)
            obj.init(cParams)
        end

        function read(obj)
            obj.openFile();
            obj.scanFile();
            obj.readDimensionAndNnode();
            obj.readCoordinates();
            obj.readConnec();
            obj.closeFile();
        end

        function d = getDataBase(obj)
            d = obj.dataBase;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.filePath = cParams.filePath;
            obj.lineNumber0 = 7;
        end

        function scanFile(obj)
            t = textscan(obj.fid,'%s','delimiter','\n', 'headerlines',0);
            obj.tScan = t{1};
        end

        function readDimensionAndNnode(obj)
            lineDimension = 5;
            tLine = obj.readLines(lineDimension);
            tLinesS = split(tLine{1});
            tLinesD = str2double(tLinesS);
            obj.dimension = tLinesD(4);
            obj.nnode = tLinesD(8);
        end

        function readCoordinates(obj)
            iline = 0;
            line0 = obj.lineNumber0;
            isCoord = true;
            xNodes = zeros(1,obj.dimension);
            while isCoord
                iline = iline + 1;
                tLine = obj.readLines(iline + line0);
                if strcmp(tLine{1},'end coordinates ')
                    isCoord = false;
                else
                    isCoord = true;
                    tLinesS = split(tLine{1});
                    tLinesD = str2double(tLinesS);
                    xNodes(iline,:) = tLinesD(2:end-1,1);
                end
            end
            obj.dataBase.coord = xNodes;
        end

        function readConnec(obj)
            line0 = obj.lineNumber0 + size(obj.dataBase.coord,1) + 3;
            iline = 0;
            isCoord = true;
            conn = zeros(1,obj.nnode);
            while isCoord
                iline = iline + 1;
                tLine = obj.readLines(iline + line0);
                if strcmp(tLine{1},'end elements')
                    isCoord = false;
                else
                    isCoord = true;
                    tLinesS = split(tLine{1});
                    tLinesD = str2double(tLinesS);
                    conn(iline,:) = tLinesD(2:obj.nnode+1,1);
                end
            end
            obj.dataBase.connec = conn;
        end

        function tLine = readLines(obj,lineNumber,nLines)
            tLine = obj.tScan(lineNumber,:);
        end

    end

end