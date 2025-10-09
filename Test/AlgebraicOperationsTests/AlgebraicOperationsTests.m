classdef AlgebraicOperationsTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
       % multiplication = {'ScalarScalar','ScalarMatrix','MatrixMatrix'}
                multiplication = {'ScalarScalar'}

        dotProduct = {'VectorVector','MatrixVector'}
        doubleDotProduct = {'MatrixMatrix','TensorMatrix'}
    end

    methods (Test, TestTags = {'Algebra'})

        function testMultiplication(testCase,multiplication)
            filename = ['testMultiplication',multiplication];
            load(filename,'input')
            [f1,f2] = testCase.createFunctions(input);
            x = f1.*f2;
            xNew = x.evaluate([0;0]);
            load(filename,'xRef');
            err = pagenorm(xNew-xRef)./pagenorm(xRef);
            tol      = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

        % function testDotProduct(testCase,dotProduct)
        %     filename = ['testDP',dotProduct];
        %     load(filename,'input')
        %     [f1,f2] = testCase.createFunctions(input);
        %     x = DP(f1,f2);
        %     xNew = x.evaluate([0;0]);
        %     load(filename,'xRef');
        %     err = pagenorm(xNew-xRef)./pagenorm(xRef); 
        %     tol      = 1e-6;
        %     testCase.verifyLessThanOrEqual(err, tol)
        % end
        % 
        % function testDoubleDotProduct(testCase,doubleDotProduct)
        %     filename = ['testDDP',doubleDotProduct];
        %     load(filename,'input')
        %     [f1,f2] = testCase.createFunctions(input);
        %     x = DDP(f1,f2);
        %     xNew = x.evaluate([0;0]);
        %     load(filename,'xRef');
        %     err = pagenorm(xNew-xRef)./pagenorm(xRef);
        %     tol      = 1e-6;
        %     testCase.verifyLessThanOrEqual(err, tol)
        % end
        
    end

    methods (Static, Access = private)
        
        function [f1,f2] = createFunctions(s)
            f1 = LagrangianFunction.create(s.mesh,s.f1.ndimf,s.f1.order);
            f1.setFValues(s.f1.values);
            f2 = LagrangianFunction.create(s.mesh,s.f2.ndimf,s.f2.order);
            f2.setFValues(s.f2.values);
        end

    end

end