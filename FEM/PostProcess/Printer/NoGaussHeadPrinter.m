classdef NoGaussHeadPrinter < HeadPrinter
    
    
    methods (Access = public)
        
        function print(obj,fileID)
            obj.fileID = fileID;
            obj.printInitialLine();
            obj.printFemMatOoHeader();
        end
        
        
    end
    
    
    
    
end