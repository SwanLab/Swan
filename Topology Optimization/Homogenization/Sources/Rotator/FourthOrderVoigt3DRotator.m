classdef FourthOrderVoigt3DRotator < FourthOrderVoigtRotator
    
    
    methods (Access = public)
        
        function obj = FourthOrderVoigt3DRotator(angle,dir)
            obj.compute(angle,dir)
        end
        
    end
    
    methods (Access = protected)
               
        function rotator = createStressRotator(a,d)
            rotator = StressVoigtRotator(a,d);
        end
        
    end

end

