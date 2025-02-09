classdef testMakeAnisotorpicTensorPlaneStressSymbolically < handle

    properties (Access = public)
        tol = 1e-10;
    end

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

        function error = computeError(obj)
            cSym = double(obj.ChSym.getValue);
            cNum = obj.ChNum.getValue;
            error = norm(cSym(:) - cNum(:));
        end

    end
    
    methods (Access = private)
        
        function createFourthOrderTensorInVoigt(obj)
            t = Stiffness3DTensor();
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

end