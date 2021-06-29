classdef ProductComputerForTensor < ProductComputer
    
    properties
    end
    
    methods (Access = public)
        
        function obj = ProductComputerForTensor(C,e)
            obj.generate(C,e)
        end
    end
    
    methods (Access = protected)
        
        function computeProduct(obj)
            C = obj.fourthOrder.getValue();
            e = obj.secondOrder.getValue();
            d = obj.secondOrder.getDimension();
            s = zeros(size(e));
            for i = 1:d
                for j = 1:d
                    for k = 1:d
                        for l = 1:d
                            cijkl = C(i,j,k,l);
                            ekl = e(k,l);
                            s(i,j) = s(i,j) + cijkl*ekl;
                        end
                    end
                end
            end
            obj.secondOrderOut.setValue(s)
        end
        
    end
    
    
end



