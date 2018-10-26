classdef testSymmetrizeFourthOrderTensor < test
    
    properties (Access = protected)
           Cani
           Csym
           En
    end
    
    methods
        
        function obj = testSymmetrizeFourthOrderTensor()
            obj.computeAnisotropicFourthOrderTensor();
            obj.computeAplicationWithRandomSymetricTensor()
        end
        
       
        
        function computeAnisotropicFourthOrderTensor(obj)
            obj.Cani = FourthOrderTensor();
            obj.Cani.createRandomTensor();
        end
        
        function computeAplicationWithRandomSymetricTensor(obj)
            Strain = StrainTensor();
            txi = Strain.tensor;
            A = obj.Cani.tensor;
            obj.En = zeros(5,1);
             for i = 1:size(A,1)
                for j = 1:size(A,2)
                    for k = 1:size(A,3)
                        for l = 1:size(A,4)
                           obj.En(1) = obj.En(1) + txi(i,j)*A(i,j,k,l)*txi(k,l);
                           obj.En(2) = obj.En(2) + txi(j,i)*A(i,j,k,l)*txi(k,l);
                           obj.En(3) = obj.En(3) + txi(i,j)*A(i,j,k,l)*txi(l,k);
                           obj.En(4) = obj.En(4) + txi(k,l)*A(i,j,k,l)*txi(i,j);
                           obj.En(5) = obj.En(5) + txi(i,k)*A(i,j,k,l)*txi(j,l);
                        end
                    end
                end
             end            
        end
        
        end
    
    methods (Access = protected)    
        function hasPassed = hasPassed(obj)
            meanEn = mean(obj.En(1:4));
            hasPassed = norm(obj.En(1:4) - meanEn) < 1e-12 && abs(meanEn - obj.En(5)) > 1e-12 ;
        end
    end
end

