classdef MatrixVectorizedInverter_3x3 < MatrixVectorizedInverter_Interface
    
    methods (Access = public)
        
        function B = computeInverse(obj,A)
            d = obj.computeAdjointsTranspose(A);
            
            det = obj.computeDeterminant(A);
            
            B = zeros(size(A));
            for i = 1:3
                for j = 1:3
                    B(i,j,:) = squeeze(d(i,j,:))./det;
                end
            end
        end
        
        function det = computeDeterminant(obj,A)
            d = obj.computeAdjointsTranspose(A);
            
            det = zeros(size(A,3),1);
            for i = 1:3
                for k = 1:3
                    det = det + squeeze(A(i,k,:).*d(k,i,:));
                end
            end
            det = det/3;
        end
        
    end
    
    methods (Access = private, Static)
        
        function d = computeAdjointsTranspose(A)
            d = zeros(size(A));
            
            d(1,1,:) = A(2,2,:).*A(3,3,:)-A(2,3,:).*A(3,2,:);
            d(1,2,:) = A(3,2,:).*A(1,3,:)-A(3,3,:).*A(1,2,:);
            d(1,3,:) = A(1,2,:).*A(2,3,:)-A(1,3,:).*A(2,2,:);
            
            d(2,1,:) = A(2,3,:).*A(3,1,:)-A(2,1,:).*A(3,3,:);
            d(2,2,:) = A(3,3,:).*A(1,1,:)-A(3,1,:).*A(1,3,:);
            d(2,3,:) = A(1,3,:).*A(2,1,:)-A(1,1,:).*A(2,3,:);
            
            d(3,1,:) = A(2,1,:).*A(3,2,:)-A(2,2,:).*A(3,1,:);
            d(3,2,:) = A(3,1,:).*A(1,2,:)-A(3,2,:).*A(1,1,:);
            d(3,3,:) = A(1,1,:).*A(2,2,:)-A(1,2,:).*A(2,1,:);
        end
        
    end
    
end