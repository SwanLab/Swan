classdef testPrintingFreeFemFile < testNotShowingError
    
    properties (Access = private)
        linesRead
        fileName
        readingFilePath        
        printingFilePath
    end
    
    methods (Access = public)
        
        function obj = testPrintingFreeFemFile()
            obj.init();            
            obj.readFile();
            obj.printFile();
        end
        
    end
    
    methods (Access = protected)
        
        function hasPassed = hasPassed(obj)
            rF = fullfile('Preprocess','FileReaders','FreeFem++',obj.readingFilePath);
            pF = fullfile(obj.printingFilePath,[obj.fileName,'.edp']);
            hasChanged = FileComparator().areFilesDifferent(rF,pF);
            hasPassed = ~hasChanged;            
        end        
        
    end    
    
    methods (Access = private)
        
        function init(obj)
            obj.fileName = 'SmoothNonSmoothRectangle';
            obj.readingFilePath = [obj.fileName,'Model','.edp'];
            obj.printingFilePath = fullfile('Input',obj.fileName);
        end
                
        function readFile(obj)
            fR = FreeFemFileReader(obj.readingFilePath);
            fR.read();
            obj.linesRead = fR.getDataBase();
        end
        
        function printFile(obj)
            d.fileName = obj.fileName;
            d.printingDir = obj.printingFilePath;  
            d.linesToPrint = obj.linesRead;
            d.type = 'Identical';
            fp = FreeFemFilePrinter.create(d);           
            fp.print();
        end
        
    end
    
end