classdef testSymmetryForIAniTensorInVoigt < test
    
    properties (Access = protected)
        SymVoigthTensor
        Ch
    end
    
    methods
        
        function obj = testSymmetryForIAniTensorInVoigt()
            obj.computeFourthOrderTensor();
            obj.computeVoigtRepresentation();
            obj.computeSymetricVoigthTensor();
        end
        
        function computeFourthOrderTensor(obj)
            obj.Ch = FourthOrderTensor();
            obj.Ch.createRandomTensor();
        end
        
        function computeVoigtRepresentation(obj)
            obj.Ch.computeTensorVoigt();
        end
        
        
        function computeSymetricVoigthTensor(obj)
            obj.SymVoigthTensor = 0.5*(obj.Ch.tensorVoigt + obj.Ch.tensorVoigt');
        end
        
    end
    
    methods (Access = protected)
        
        function hasPassed = hasPassed(obj)
            hasPassed = norm(double(obj.SymVoigthTensor(:)) - obj.Ch.tensorVoigt(:)) < 1e-6;
        end
        
    end
end

