classdef testSymmetrizeFourthOrderTensor < handle
    
    properties
           Cani
           Csym
           fileName = 'testSymmetrizeFourthOrderTensor'
           En
    end
    
    methods
        
        function obj = testSymmetrizeFourthOrderTensor()
            obj.printHeader()
            obj.computeAnisotropicFourthOrderTensor();
            obj.computeAplicationWithRandomSymetricTensor()
            obj.checkTestPassed()
            obj.printTail()
        end
        
       
        
        function computeAnisotropicFourthOrderTensor(obj)
            obj.Cani = fourthOrderTensor();
            obj.Cani.createRandomTensor();
            obj.Cani.MakeMajorAndMinorSymmetrization()
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
        
        
        function hasPassed = hasPassed(obj)
            meanEn = mean(obj.En(1:4));
            hasPassed = norm(obj.En(1:4) - meanEn) < 1e-12 && abs(meanEn - obj.En(5)) > 1e-12 ;
        end
        
        function checkTestPassed(obj)
            if obj.hasPassed()
                cprintf('green',strcat(obj.fileName,' PASSED\n'));
            else
                cprintf('err',strcat(obj.fileName,' FAILED\n'));
            end
        end
        
        
    end
    
    methods (Static)
        
        function printHeader()
            fprintf('Running TopOpt tests...\n')
        end
        
        function printTail()
            fprintf('\nTopOpt tests completed.\n')
            fprintf('\n-------------------------------------------\n\n')
        end
        

    end
end

