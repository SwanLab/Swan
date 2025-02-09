classdef SecondOrder3DVoigtTensor < AbstractTensor ...
                                  & SecondOrderDescriptor ...
                                  & VoigtRepresentation ...
                                  & Elasticity3dDescriptor
    
    methods (Access = public)
        
        function obj = SecondOrder3DVoigtTensor()
        end
    end
    
    methods (Access = protected)
        
        function loadTensorSize(obj)
            obj.tensorSize = [6,1];
        end
    end
    
    
    
end

