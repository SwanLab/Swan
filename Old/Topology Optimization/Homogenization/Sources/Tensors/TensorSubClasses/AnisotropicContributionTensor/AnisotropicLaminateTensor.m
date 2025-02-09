classdef AnisotropicLaminateTensor < SymmetricFourthOrder3DTensor
    
    properties (Access = protected)
        direction
        materialTensor
        mu
    end
    
    methods (Access = protected)
        
        function init(obj,A,dir)
            obj.direction      = dir;
            obj.materialTensor = A.getValue();
            obj.mu             = A.getMu();
        end
        
        function t = createNullFourthOrderTensor(obj)
            t = zeros(obj.getTensorSize);
            if obj.isSymbollic()
                t = sym(t);
            end
        end
        
        function itIs = isSymbollic(obj)
            isTensorSym    = isa(obj.materialTensor,'sym');
            isDirectionSym = isa(obj.direction,'sym');
            itIs = isTensorSym || isDirectionSym;
        end

    end
    
end

