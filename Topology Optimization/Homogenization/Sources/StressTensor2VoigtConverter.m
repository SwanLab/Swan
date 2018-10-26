classdef StressTensor2VoigtConverter < SecondOrderTensor2VoigtConverter
    
    properties
    end
    
    methods (Access = public)
        
        function obj = StressTensor2VoigtConverter(tensor)
           obj.computeConversion(tensor) 
        end
        
    end
    
    methods (Access = protected)
        
        function factor = computeVoigtFactor(obj)
            factor = 1;
        end
    end
    
end

