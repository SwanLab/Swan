classdef SecondOrderVoigt2TensorConverterPS < Voigt2TensorConverter ...
                                            & Voigt2TensorPSrepresentation ...
                                            & SecondOrderVoigt2TensorDescriptor
    
    properties (Access = protected)
        voigtIndex
    end
      
    methods (Access = protected)
             
        function obtainTensorSize(obj)
            obj.tensorSize = [obj.dim,obj.dim];
        end
        
    end
    

    
end