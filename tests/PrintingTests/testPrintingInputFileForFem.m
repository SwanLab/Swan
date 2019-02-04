classdef testPrintingInputFileForFem < testNotShowingError
    
    properties (Access = private)
        fileName 
        gmsFile 
        resultsDir
        fullFileName
    end
    
    methods (Access = public)
        
        function obj = testPrintingInputFileForFem()
            obj.init()
            g  = obj.gmsFile;
            oD = obj.resultsDir;
            oF = obj.fullFileName;
            c = GmsFile2FemMatOoFileConverter(g,oD,oF);
            c.convert();
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
            obj.fileName = 'InputFileForFem';
            obj.gmsFile = 'testReadingGmsh.msh';            
            obj.resultsDir   = fullfile('Output',obj.fileName);
            obj.fullFileName = fullfile(obj.resultsDir,[obj.fileName,'.m']);
        end

    end
    
    
end
