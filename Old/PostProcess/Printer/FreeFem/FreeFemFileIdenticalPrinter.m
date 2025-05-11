classdef FreeFemFileIdenticalPrinter < FreeFemFilePrinter
        
    methods (Access = protected)
                       
       function printingLines(obj)
            nlines = numel(obj.linesToPrint);
            for iline = 1:nlines
                line = obj.linesToPrint{iline};               
                fprintf(obj.fileID,'%s\r\n',line);
            end
        end
    
    end
    
end