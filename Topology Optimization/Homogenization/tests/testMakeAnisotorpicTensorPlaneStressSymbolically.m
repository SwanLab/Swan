classdef testMakeAnisotorpicTensorPlaneStressSymbolically < test
    
    properties (Access = protected)
        ChVoigt
        ChSym
        ChNum
    end
    
    methods (Access = public)
        
        function obj = testMakeAnisotorpicTensorPlaneStressSymbolically()
           obj.createFourthOrderTensorInVoigt()
           obj.createPlaneStressTensorSymbolically()
           obj.createNumericFourthOrderTensor();
        end
        
    end
    
    methods (Access = private)
        
        function createFourthOrderTensorInVoigt(obj)
            t = SymmetricFourthOrder3DTensor();
            t.createRandomTensor();
            obj.ChVoigt = Tensor2VoigtConverter.convert(t);
        end
        
        function createPlaneStressTensorSymbolically(obj)
            t = obj.ChVoigt.getValue();
            tPS = PST4VoigtFourthOrderTensorSymbolically(t);            
            obj.ChSym = tPS;
        end
        
        function createNumericFourthOrderTensor(obj)
           t = obj.ChVoigt;
           obj.ChNum = PlaneStressTransformer.transform(t); 
        end
    end
    
    
    methods (Access = protected)
        function hasPassed = hasPassed(obj)
            cSym = obj.ChSym.getValue;
            cNum = obj.ChNum.getValue;
            hasPassed = norm(cSym(:) - cNum(:)) < 1e-13;
        end
    end
end


