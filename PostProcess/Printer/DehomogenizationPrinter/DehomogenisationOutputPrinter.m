classdef DehomogenisationOutputPrinter < FilePrinter
    
    properties (Access = private)
        values
    end
    
    methods (Access = public)
        
        function obj = DehomogenisationOutputPrinter(cParams)
            obj.init(cParams)            
        end
        
        function print(obj)
            obj.openFile();
            obj.printLines();
            obj.closeFile();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName = cParams.fileName;
            obj.values   = cParams.values;
        end
        
        function printLines(obj)
           formatV = repmat('%4.12f ',1,size(obj.values,2));
           format = [formatV,'\n'];
           fprintf(obj.fileID,format,obj.values');
        end
        
    end
    
end