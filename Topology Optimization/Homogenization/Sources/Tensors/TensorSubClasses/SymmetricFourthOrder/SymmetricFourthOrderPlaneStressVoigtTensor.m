classdef SymmetricFourthOrderPlaneStressVoigtTensor <  AbstractTensor ...
                                              & FourthOrderDescriptor ...
                                              & VoigtRepresentation ...
                                              & ElasticityPlaneStressDescriptor
    
    methods (Access = public)
        
        function obj = SymmetricFourthOrderPlaneStressVoigtTensor()
        end
        
        function d = getVoigtDimension(obj)
            d = obj.dimVoigt;
        end
    end
    
    methods (Access = protected)
                
        function loadTensorSize(obj)
            obj.tensorSize = [3,3];
        end
    end
    
end

