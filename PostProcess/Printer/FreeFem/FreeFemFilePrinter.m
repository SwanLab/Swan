classdef FreeFemFilePrinter < FilePrinter
    
    properties (Access = private)
        printingDir
    end
    
    properties (Access = protected)
        linesToPrint        
    end
    
    methods (Access = public, Static)
        
        function p = create(inputData)
            p = FreeFemFilePrinterFactory.create(inputData);
            p.init(inputData);
            p.createOutPutFolder();            
        end
        
    end
    
    methods (Access = public)
        
        function print(obj)
            obj.openFile();
            obj.printingLines();
            obj.closeFile();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,d)
            obj.linesToPrint = d.linesToPrint;
            obj.printingDir  = d.printingDir;
            obj.fileName     = fullfile(obj.printingDir,[d.fileName,'.edp']);            
        end        
                
        function createOutPutFolder(obj)
            dir = obj.printingDir;
            if ~exist(dir,'dir')
                mkdir(dir)
                addpath(dir)
            end
        end
        
    end
    
    methods (Access = protected, Abstract)
        printingLines(obj)        
    end
end