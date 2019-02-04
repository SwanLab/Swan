classdef GmsFile2FemMatOoFileConverter < handle
    
    properties (Access = private)
        readData         
        gmsFile
        outPutDir
        outPutFileName
    end
    
    methods (Access = public)
        
        function obj = GmsFile2FemMatOoFileConverter(gF,oD,oF)
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
            rD.connec = reader.connec;
            rD.coord  = reader.coord;
            rD.isElemInThisSet = reader.isElemInThisSet;
            rD.masterSlave = reader.masterSlave;
            rD.corners = reader.corners;
            obj.readData = rD;
        end
        
        function printInputFemFile(obj)
            data = obj.readData;
            data.resultsDir = obj.outPutDir;
            data.fileName   = obj.outPutFileName;
            InputFemFilePrinter(data);
        end        
        
    end
    
    
end