classdef MatrixVectorizedInverter_2x2 < MatrixVectorizedInverter_Interface
    
    methods (Access = public)
        
        function B = computeInverse(obj,A)
            d(1,1,:,:) = A(2,2,:,:);
            d(1,2,:,:) = -A(1,2,:,:);
            d(2,1,:,:) = -A(2,1,:,:);
            d(2,2,:,:) = A(1,1,:,:);
            
            det = obj.computeDeterminant(A);
            
            B = zeros(size(A));
            for i = 1:2
                for j = 1:2
                    B(i,j,:,:) = squeeze(d(i,j,:,:))./det;
                end
            end
        end
        
        function detA = computeDeterminant(~,A)
            detA = squeeze(A(1,1,:,:).*A(2,2,:,:)-A(1,2,:,:).*A(2,1,:,:));
        end
        
    end
    
end