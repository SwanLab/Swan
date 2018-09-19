classdef testSymmetrizeIsotropicFourthOrderTensor < handle
    
    properties
           Ciso
           Csym
           fileName = 'testSymmetrizeIsotropicFourthOrderTensor'
    end
    
    methods
        
        function obj = testSymmetrizeIsotropicFourthOrderTensor()
            obj.printHeader()
            obj.computeIsotropicFourthOrderTensor();
            obj.creatreSymmetricFourthOrderTensor()
            obj.checkTestPassed()
            obj.printTail()
        end
        
       
        
        function computeIsotropicFourthOrderTensor(obj)
            E = 1; nu = 1/3;
            obj.Ciso = IsotropicConstitutiveTensor3D(E,nu);
        end
        
        function creatreSymmetricFourthOrderTensor(obj)
            obj.Csym = obj.Ciso;
            obj.Csym.MakeMajorAndMinorSymmetrization();            
        end
        
        
        function hasPassed = hasPassed(obj)
            hasPassed = norm(obj.Csym.tensor(:) - obj.Ciso.tensor(:)) < 1e-6;
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

