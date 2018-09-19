classdef testSymmetryForIAniTensorInVoigt < handle
    
    properties
        fileName = 'testSymmetryForIAniTensorInVoigt';
        SymVoigthTensor
        Ch
    end
    
    methods
        
        function obj = testSymmetryForIAniTensorInVoigt()
            obj.printHeader()
            obj.computeFourthOrderTensor();
            obj.computeVoigtRepresentation();
            obj.computeSymetricVoigthTensor();
            obj.checkTestPassed()
            obj.printTail()
        end
        
        function computeFourthOrderTensor(obj)
            obj.Ch = fourthOrderTensor();
            obj.Ch.createRandomTensor();
            obj.Ch.MakeMajorAndMinorSymmetrization()
        
        end
        
        function computeVoigtRepresentation(obj)
                obj.Ch.tensorVoigt = obj.Ch.RespresentTensorinVoigt(obj.Ch.tensor);
        end
        
        
        function computeSymetricVoigthTensor(obj)
                obj.SymVoigthTensor = 0.5*(obj.Ch.tensorVoigt + obj.Ch.tensorVoigt');
        end
        
        
        
        function hasPassed = hasPassed(obj)
            hasPassed = norm(double(obj.SymVoigthTensor(:)) - obj.Ch.tensorVoigt(:)) < 1e-6;
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

