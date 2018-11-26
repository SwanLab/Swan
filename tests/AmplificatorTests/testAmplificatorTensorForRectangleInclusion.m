classdef testAmplificatorTensorForRectangleInclusion < ...
         testAmplificatorTensorForInclusions
    
    
    properties (Access = protected)
        testName = 'AmplificatorTensorForRectangleInclusion';
    end
    
    methods (Access = protected)
        
        function h = createHomogenizer(obj,iter,m1)
            f     = obj.fileName;
            m2    = 0.2;
            print = obj.printTopology;
            h     = NumericalRectangleHomogenizer(f,print,m1,m2,iter);
        end
        
    end    

end

