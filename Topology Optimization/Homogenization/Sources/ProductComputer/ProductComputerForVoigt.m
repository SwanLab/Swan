classdef ProductComputerForVoigt < ProductComputer
    
    properties
    end
    
    methods (Access = public)
        
        function obj = ProductComputerForVoigt(C,e)
            obj.generate(C,e)
        end
    end
    
    methods (Access = protected)
        
        function computeProduct(obj)
            C = obj.fourthOrder.getValue();
            e = obj.secondOrder.getValue();
            d = obj.secondOrder.getVoigtDimension();
            s = zeros(size(e));
            for i = 1:d
                for j = 1:d
                    cij = C(i,j);
                    ej = e(j);
                    s(i) = s(i) + cij*ej;
                end
            end
            obj.secondOrderOut.setValue(s)
        end
        
    end
    
    
end



