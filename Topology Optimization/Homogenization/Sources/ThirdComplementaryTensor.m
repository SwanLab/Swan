classdef ThirdComplementaryTensor < FourthOrderTensor
    
    properties
        direction
    end
    
    methods
        
        function obj = ThirdComplementaryTensor(A,direction)
            obj.direction = direction;
            obj.generate(A)
        end
        
        function  generate(obj,A)
            obj.tensor = obj.compute(A);
        end
        
        function F2ten = compute(obj,IsoTensor)
            T = IsoTensor.getValue();
            mu = IsoTensor.getMu();
            lambda = IsoTensor.getLambda();
            
            dim = size(T,1);
            F2ten =  sym(zeros([dim dim dim dim]));
            
            for i = 1:dim
                for j = 1:dim
                    
                    T1 = sym(zeros([dim dim]));
                    for p = 1:dim
                        for q = 1:dim
                            T1(i,j) = T1(i,j) + T(p,q,i,j)*obj.direction(p)*obj.direction(q);
                        end
                    end
                    
                    F1v = sym(zeros([dim 1]));
                    for p = 1:dim
                        for q = 1:dim
                            F1v(q) = F1v(q) +  T(i,j,p,q)*obj.direction(p);
                        end
                    end
                    
                    for k = 1:dim
                        for l = 1:dim
                            T2 = sym(zeros([dim dim]));
                            for p = 1:dim
                                for q = 1:dim
                                    T2(k,l) = T2(k,l) + T(p,q,k,l)*obj.direction(p)*obj.direction(q);
                                end
                            end
                            F2 = (mu+lambda)/(mu*(2*mu+lambda))*T1(i,j)*T2(k,l);
                            F2ten(i,j,k,l) = F2;

                        end
                    end
                end
            end
            
            
        end
        

    end
    
end

