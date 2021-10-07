classdef TestCarpet < handle & matlab.unittest.TestCase
    % Una suite agrupa tests  els fa correr de manera relativament
    % senzilla. Per exemple, aqui sota, hi ha un testcase construit amb
    % diferents parametres per generar tests
    properties (TestParameter)
        type = {'single','double','uint16'};
        level = struct('small',2,'medium',4,'large',6);
        side = struct('small',9,'medium',81,'large',729);
    end

    methods (Test) % En total n'hi ha 16
        function testRemainPixels(testCase,level) % 3 tests (3 levels)
            expPixelCount = 8^level;
            actPixels = find(sierpinski(level));
            testCase.verifyNumElements(actPixels,expPixelCount)
        end
        
        function testClass(testCase,type,level) % 9 tests
            testCase.verifyClass(sierpinski(level,type),type)
        end
        
        function testDefaultL1Output(testCase) % 1 test
            exp = single([1 1 1; 1 0 1; 1 1 1]);
            testCase.verifyEqual(sierpinski(1),exp)
        end
    end
    
    methods (Test, ParameterCombination = 'sequential')
        function testNumel(testCase,level,side) % 3 tests
            import matlab.unittest.constraints.HasElementCount
            testCase.verifyThat(sierpinski(level), ...
                HasElementCount(side^2))
        end
    end 
end