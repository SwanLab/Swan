classdef testSymmetryForIAniTensorInVoigt < testShowingError
    
    properties (Access = protected)
        ChVoigtSym
        Ch
        ChVoigt
        
        tol = 1e-12;
    end
    
    methods
        
        function obj = testSymmetryForIAniTensorInVoigt()
            obj.computeFourthOrderTensor();
            obj.computeVoigtRepresentation();
            obj.computeSymetricVoigthTensor();
        end
        
        function computeFourthOrderTensor(obj)
            obj.Ch = Stiffness3DTensor();
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
        
        function computeError(obj)
            c    = obj.ChVoigt.getValue();
            cSym = obj.ChVoigtSym;            
            obj.error = norm(c(:) - cSym(:));
        end
        
    end
end
