classdef testSymmetryForIAniTensorInVoigt < handle
    
    properties (Access = protected)
        ChVoigtSym
        Ch
        ChVoigt
    end

    properties (Access = public)
        tol = 1e-12;
    end
    
    methods (Access = public)
        
        function obj = testSymmetryForIAniTensorInVoigt()
            obj.computeFourthOrderTensor();
            obj.computeVoigtRepresentation();
            obj.computeSymetricVoigthTensor();
        end
        
        function error = computeError(obj)
            c    = obj.ChVoigt.getValue();
            cSym = obj.ChVoigtSym;
            error = norm(c(:) - cSym(:));
        end

    end

    methods (Access = private)
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

end