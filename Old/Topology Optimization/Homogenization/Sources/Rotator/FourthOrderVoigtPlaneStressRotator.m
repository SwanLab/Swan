classdef FourthOrderVoigtPlaneStressRotator < FourthOrderVoigtRotator

    methods (Access = public)
        
        function obj = FourthOrderVoigtPlaneStressRotator(angle,dir)
            obj.compute(angle,dir)
        end
        
    end
    
    methods (Access = protected,Static)
        
        function r = createStressRotator(a,d)
            r = StressVoigtPlaneStressRotator(a,d);
        end
       
    end
    
end

