classdef GiDPrinter < FilePrinter
    
    properties (Access = protected, Abstract)
       ext 
    end
    
    properties (Access = protected)
        iter
        ndim
        etype
        outFileName
        resultsDir
    end
    
    
    methods (Access = protected)
        
       function createFileName(obj,iS)
            fileNameWithNoPath = strcat(obj.outFileName,num2str(iS),'.flavia.',obj.ext);
            obj.fileName = fullfile(obj.resultsDir,fileNameWithNoPath);
       end
        
    end
    
    
end