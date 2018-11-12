classdef ThirdAnisotropicTensor < AnisotropicLaminateTensor
    
    properties (Access = private)
        lambda
    end
    
    methods (Access = public)
        
        function obj = ThirdAnisotropicTensor(A,dir)
            obj.init(A,dir)
            obj.compute()
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,A,dir)
            obj.init@AnisotropicLaminateTensor(A,dir)
            obj.lambda = A.getLambda();
        end
        
    end
    
    methods (Access = private)
        
        function compute(obj)
            t = obj.materialTensor;
            muV = obj.mu;
            lam = obj.lambda;
            factor = (muV+lam)/(muV*(2*muV+lam));
            dir = obj.direction;            
            d = obj.getDimension;
            A = obj.createNullFourthOrderTensor();  
            for i = 1:d
                for j = 1:d
                    T1 = obj.computeSecondOrderTensor(t,dir,i,j);
                    for k = 1:d
                        for l = 1:d
                            T2 = obj.computeSecondOrderTensor(t,dir,k,l);                             
                            aijkl = factor*T1(i,j)*T2(k,l);
                            A(i,j,k,l) = aijkl;
                        end
                    end
                end
            end            
            obj.setValue(A)            
        end
        
        function a = computeSecondOrderTensor(obj,t,dir,i,j)
            d = obj.getDimension;
            a = obj.createNullSecondOrderTensor();
            for p = 1:d
                for q = 1:d
                    a(i,j) = a(i,j) + t(p,q,i,j)*dir(p)*dir(q);
                end
            end
        end
        
        function t = createNullSecondOrderTensor(obj)
            dim = obj.getDimension;
            t = zeros(dim,dim);
            if obj.isSymbollic()
                t = sym(t);
            end
        end
        
    end
    
end

