classdef NoGaussHeadPrinter < HeadPrinter
    
    methods (Access = public)
        
        function obj = NoGaussHeadPrinter(dh)
           obj.fileID = dh.fileID;
        end
        
        function print(obj)
            obj.printInitialLine();
            obj.printFemMatOoHeader();
        end
        
    end
    
end