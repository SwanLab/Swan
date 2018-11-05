classdef testSymmetryForIAniTensorInVoigt < test
    
    properties (Access = protected)
        ChVoigtSym
        Ch
        ChVoigt
    end
    
    methods
        
        function obj = testSymmetryForIAniTensorInVoigt()
            obj.computeFourthOrderTensor();
            obj.computeVoigtRepresentation();
            obj.computeSymetricVoigthTensor();
        end
        
        function computeFourthOrderTensor(obj)
            obj.Ch = SymmetricFourthOrder3DTensor();
            obj.Ch.createRandomTensor();
        end
        
        function computeVoigtRepresentation(obj)
            obj.ChVoigt = Tensor2VoigtConverter.convert(obj.Ch);
        end
        
        
        function computeSymetricVoigthTensor(obj)
            t = obj.ChVoigt.getValue();
            obj.ChVoigtSym = 0.5*(t + t');
        end
        
    end
    
    methods (Access = protected)
        
        function hasPassed = hasPassed(obj)
            c    = obj.ChVoigt.getValue();
            cSym = obj.ChVoigtSym;            
            hasPassed = norm(c(:) - cSym(:)) < 1e-12;
        end
        
    end
end

