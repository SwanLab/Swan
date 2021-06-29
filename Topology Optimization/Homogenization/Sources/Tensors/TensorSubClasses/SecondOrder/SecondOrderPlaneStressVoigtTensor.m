classdef SecondOrderPlaneStressVoigtTensor < AbstractTensor ...
                                        & SecondOrderDescriptor ...
                                        & VoigtRepresentation ...
                                        & ElasticityPlaneStressDescriptor
    
    methods (Access = public)
        
        function obj = SecondOrderPlaneStressVoigtTensor()
        end
    end
            
    methods (Access = protected)
        
        function loadTensorSize(obj)
            obj.tensorSize = [3,1];
        end
    end
    
end

