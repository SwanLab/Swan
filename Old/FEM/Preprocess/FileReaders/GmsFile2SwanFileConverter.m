classdef GmsFile2SwanFileConverter < handle
    
    properties (Access = private)
        readData
        gmsFile
        outPutDir
        outPutFileName
    end
    
    methods (Access = public)
        
        function obj = GmsFile2SwanFileConverter(gF,oD,oF)
            obj.init(gF,oD,oF)
        end
        
        function convert(obj)
            obj.readGmsFile();
            obj.printInputFemFile();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,gmsFile,outPutDir,outPutFileName)
            obj.gmsFile = gmsFile;
            obj.outPutDir = outPutDir;
            obj.outPutFileName = outPutFileName;
        end
        
        function readGmsFile(obj)
            reader = GmsReader(obj.gmsFile);
            reader.read();
            rDB = reader.getDataBase();
            rD.connec = rDB.connec;
            rD.coord  = rDB.coord;
            rD.isElemInThisSet = rDB.isElemInThisSet;
            rD.masterSlave = rDB.masterSlave;
            rD.corners = rDB.corners;
            obj.readData = rD;
        end
        
        function printInputFemFile(obj)
            data = obj.readData;
            data.resultsDir = obj.outPutDir;
            data.fileName   = obj.outPutFileName;
            fp = InputFemFilePrinter(data);
            fp.print();
        end
        
    end
    
end