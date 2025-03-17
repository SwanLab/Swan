classdef SymmetricFourthOrder3DVoigtTensor <  AbstractTensor ...
                                              & FourthOrderDescriptor ...
                                              & VoigtRepresentation ...
                                              & Elasticity3dDescriptor
    
                                          
    methods (Access = public)
        
        function obj = SymmetricFourthOrder3DVoigtTensor()
        end
        
        function createRandomTensor(obj)
            obj.createRandomTensor@AbstractTensor();
            obj.makeSymmetrization();
        end
    end
    
    methods (Access = protected)
        
        function makeSymmetrization(obj)
            t = obj.getValue;
            ts = 0.5*(t + t');
            obj.setValue(ts);
        end
        
        function loadTensorSize(obj)
            obj.tensorSize = [6,6];
        end
    end
    
end

