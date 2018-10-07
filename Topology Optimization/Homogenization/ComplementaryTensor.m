classdef ComplementaryTensor < fourthOrderTensor
    
    properties
        IsotropicConstitutiveTensor
        direction
    end
    
    methods
        
        function obj = ComplementaryTensor(A,direction)
            obj.IsotropicConstitutiveTensor = A; 
            obj.direction = direction;
            obj.generate()
            obj.tensorVoigt = obj.RespresentTensorinVoigt(obj.tensor);
        end
        
       function  generate(obj)
           dim = 3;
           IsoTensor = obj.IsotropicConstitutiveTensor.tensor;
           mu = obj.IsotropicConstitutiveTensor.mu;
           lambda = obj.IsotropicConstitutiveTensor.kappa; 
           
           %obj.tensor = zeros(dim,dim,dim,dim);
           obj.tensor = sym('tensor',[dim dim dim dim]);
           
               for i = 1:dim
                    for j = 1:dim
                        
                        T1 = sym(zeros([dim dim])); 
                        for p = 1:dim
                            for q = 1:dim
                                T1(i,j) = T1(i,j) + IsoTensor(p,q,i,j)*obj.direction(p)*obj.direction(q);
                            end
                        end
                        
                        F1v = sym(zeros([dim 1]));
                        for p = 1:dim
                            for q = 1:dim
                                F1v(q) = F1v(q) +  IsoTensor(i,j,p,q)*obj.direction(p);
                            end
                        end
                       
                        for k = 1:dim
                            for l = 1:dim
                                
                                F0 = IsoTensor(i,j,k,l);
                                
                                %F1v = zeros(dim,1);
                                
                                F2v = sym(zeros([dim 1]));
                                for p = 1:dim
                                    for q = 1:dim
                                        F2v(q) = F2v(q) +  IsoTensor(k,l,p,q)*obj.direction(p);
                                    end
                                end
                                
                                
                                F1t = sum(F1v(:).*F2v(:));
                                F1 = -1/mu*F1t;
                                
                                T2 = sym(zeros([dim dim])); 
                                for p = 1:dim
                                    for q = 1:dim
                                        T2(k,l) = T2(k,l) + IsoTensor(p,q,k,l)*obj.direction(p)*obj.direction(q);
                                    end
                                end
                                F2 = (mu+lambda)/(mu*(2*mu+lambda))*T1(i,j)*T2(k,l);
                                
                                obj.tensor(i,j,k,l) = F0 + F1 + F2;
                            end
                        end
                    end
                end
        end
        
        
    end
    
end

