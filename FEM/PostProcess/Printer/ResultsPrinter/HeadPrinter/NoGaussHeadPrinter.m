classdef NoGaussHeadPrinter < HeadPrinter
    
    
    methods (Access = public)
        
        function print(obj,hD)
            obj.fileID = hD.fileID;
            obj.printInitialLine();
            obj.printFemMatOoHeader();
        end
        
        
    end
    
    
    
    
end