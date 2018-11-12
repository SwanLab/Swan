classdef SecondAnisotropicTensor < AnisotropicLaminateTensor
    

    
    methods (Access = public)
        
        function obj = SecondAnisotropicTensor(A,dir)
            obj.init(A,dir)
            obj.compute()
        end
        
    end
    
    methods (Access = private)
        
        function  compute(obj)
            t = obj.materialTensor;
            dir = obj.direction;
            d = obj.getDimension;
            A = obj.createNullFourthOrderTensor();
            for i = 1:d
                for j = 1:d
                    v1 = obj.computeVector(t,dir,i,j);
                    for k = 1:d
                        for l = 1:d                            
                            v2 = obj.computeVector(t,dir,k,l);                        
                            v1v2 = sum(v1(:).*v2(:));
                            aijkl = -1/obj.mu*v1v2;
                            A(i,j,k,l) = aijkl;                            
                        end
                    end
                end
            end            
            obj.setValue(A)            
        end
        
        function v1 = computeVector(obj,t,dir,i,j)
            d = obj.getDimension;
            v1 = obj.createNullVector();
            for p = 1:d
                for q = 1:d
                    v1(q) = v1(q) +  t(i,j,p,q)*dir(p);
                end
            end
        end        
        
        function v = createNullVector(obj)
            v = zeros(obj.getDimension,1);
            if obj.isSymbollic()
                v = sym(v);
            end
        end
        
    end
    
end

