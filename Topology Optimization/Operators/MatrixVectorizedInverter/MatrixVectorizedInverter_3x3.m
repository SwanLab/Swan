classdef MatrixVectorizedInverter_3x3 < MatrixVectorizedInverter_Interface
    
    methods (Access = public)
        
        function B = computeInverse(obj,A)
            d = obj.computeAdjointsTranspose(A);
            det = obj.computeDeterminant(A);
            B = zeros(size(A));
            for i = 1:3
                for j = 1:3
                    B(i,j,:,:) = squeeze(d(i,j,:,:))./det;
                end
            end
        end
        
        function detA = computeDeterminant(obj,A)
            d    = obj.computeAdjointsTranspose(A);
            d = permute(d,[2 1 3 4]);
            detA = squeeze(sum(A.*d,[1 2]))/3;
        end
        
    end
    
    methods (Access = private, Static)
        
        function d = computeAdjointsTranspose(A)
            d = zeros(size(A));
            
            A11 = A(1,1,:,:);
            A12 = A(1,2,:,:);
            A13 = A(1,3,:,:);
            A21 = A(2,1,:,:);
            A22 = A(2,2,:,:);
            A23 = A(2,3,:,:);
            A31 = A(3,1,:,:);
            A32 = A(3,2,:,:);
            A33 = A(3,3,:,:);
           
            d(1,1,:,:) = A22.*A33-A23.*A32;
            d(1,2,:,:) = A32.*A13-A33.*A12;
            d(1,3,:,:) = A12.*A23-A13.*A22;
            
            d(2,1,:,:) = A23.*A31-A21.*A33;
            d(2,2,:,:) = A33.*A11-A31.*A13;
            d(2,3,:,:) = A13.*A21-A11.*A23;
            
            d(3,1,:,:) = A21.*A32-A22.*A31;
            d(3,2,:,:) = A31.*A12-A32.*A11;
            d(3,3,:,:) = A11.*A22-A12.*A21;
        end
        
    end
    
end