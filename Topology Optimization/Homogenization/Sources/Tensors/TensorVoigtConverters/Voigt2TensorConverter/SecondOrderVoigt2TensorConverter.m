classdef SecondOrderVoigt2TensorConverter < Voigt2TensorConverter ...
                                            & Voigt2Tensor3Drepresentation ...
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