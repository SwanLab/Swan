classdef testPrintingInputFileForFem < testNotShowingError
    
    properties (Access = private)
        readData
        fileName = 'InputFileForFem'
        resultsDir
        fullFileName
    end
    
    methods (Access = public)
        
        function obj = testPrintingInputFileForFem()
            obj.init()
            obj.readGmsFile();
            obj.printInputFemFile();
        end
        
    end
    
    methods (Access = protected)
        
        function hasPassed = hasPassed(obj)
            sF = fullfile('tests','PrintingTests','PrintedFiles',['test',obj.fileName,'.m']);
            oF = obj.fullFileName;
            hasChanged = FileComparetor.areFilesDifferent(sF,oF);
            hasPassed = ~hasChanged;
        end        
        
    end
    
    
    methods (Access = private)
        
        function init(obj)
            obj.resultsDir = fullfile('Output',obj.fileName);
            obj.fullFileName = fullfile(obj.resultsDir,[obj.fileName,'.m']);
        end
        
        function readGmsFile(obj)
           filePath = 'tests/ReadingFilesTests/ReadingFiles/testReadingGmsh.msh';
           reader = GmsReader(filePath);
           rD.connec = reader.connec;
           rD.coord  = reader.coord;
           rD.isElemInThisSet = reader.isElemInThisSet;
           rD.masterSlave = reader.masterSlave;
           rD.corners = reader.corners;           
           obj.readData = rD;
        end
        
        function printInputFemFile(obj)
            data = obj.readData;
            data.resultsDir = obj.resultsDir;
            data.fileName   = obj.fullFileName;
            InputFemFilePrinter(data);
        end
        
    end
    
    
end
