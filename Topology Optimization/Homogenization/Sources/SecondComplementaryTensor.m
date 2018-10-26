classdef SecondComplementaryTensor < FourthOrderTensor
    
    properties
        direction
    end
    
    methods
        
        function obj = SecondComplementaryTensor(A,direction)
            obj.direction = direction;
            obj.generate(A)
        end
        
        function  generate(obj,A)
            obj.tensor = obj.compute(A);
        end
        
        function F1ten = compute(obj,Tensor)
            IsoTensor = Tensor.getValue;
            mu = Tensor.getMu;
            
            dim = size(IsoTensor,1);
            F1ten =  sym(zeros([dim dim dim dim]));
            
            for i = 1:dim
                for j = 1:dim
                    
                    F1v = sym(zeros([dim 1]));
                    for p = 1:dim
                        for q = 1:dim
                            F1v(q) = F1v(q) +  IsoTensor(i,j,p,q)*obj.direction(p);
                        end
                    end
                    
                    for k = 1:dim
                        for l = 1:dim
                            
                            F2v = sym(zeros([dim 1]));
                            for p = 1:dim
                                for q = 1:dim
                                    F2v(q) = F2v(q) +  IsoTensor(k,l,p,q)*obj.direction(p);
                                end
                            end
                            
                            F1t = sum(F1v(:).*F2v(:));
                            F1 = -1/mu*F1t;
                            F1ten(i,j,k,l) = F1;
                            
                        end
                    end
                end
            end
            
        end
        
    end
    
end

