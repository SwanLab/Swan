classdef FourthWithSecondOrderProductComputer
    
    properties
    end
    
    methods (Access = public)
        
        function obj = FourthWithSecondOrderProductComputer()
            
        end
        
    end
    
    methods (Static, Access = public)        

        function tensor = computeInTensor(fourthTensor,secondTensor)
            dim = size(secondTensor,1);
            tensor = zeros(size(secondTensor));
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            fourth = fourthTensor(i,j,k,l);
                            second = secondTensor(k,l);
                            tensor(i,j) = tensor(i,j) + fourth*second;
                        end
                    end
                end
            end
        end
        
        function tensor = computeInVoigt(fourthTensor,secondTensor)
            dim = size(secondTensor,1);
            tensor = zeros(size(secondTensor));
            for i = 1:dim
                for j = 1:dim
                    fourth = fourthTensor(i,j);
                    second = secondTensor(j);
                    tensor(i) = tensor(i) + fourth*second;
                end
            end
        end
        
    end
    
end

