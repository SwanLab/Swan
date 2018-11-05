classdef SymmetricFourthOrder3DVoigtTensor <  AbstractTensor ...
                                              & FourthOrderDescriptor ...
                                              & VoigtRepresentation ...
                                              & Elasticity3dDescriptor
    
                                          
    methods (Access = public)
        
        function obj = SymmetricFourthOrder3DVoigtTensor()
        end
    end
    
    methods (Access = protected)
        
        function loadTensorSize(obj)
            obj.tensorSize = [6,6];
        end
    end
    
end

