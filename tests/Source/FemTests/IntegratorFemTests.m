classdef IntegratorFemTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        triangle = {'test2d_triangle'}
    end

    methods (Test, TestTags = {'Triangle'})

%         function testWholeTriangle(testCase, triangle)
%             s.computerType    = 'FEM';
%             s.testName         = triangle;
%             s.variablesToStore = {'d_u'};
%             test = PrecomputedVariableTest(s);
%             err = test.computeError();
%             tol = 1e-6;
%             testCase.verifyLessThanOrEqual(err, tol)
%         end

        function testLHS(testCase, triangle)
            s.computerType    = 'FEM';
            s.testName         = triangle;
            s.variablesToStore = {'d_u'};
            test = PrecomputedVariableTest(s);
            comp = test.getComputation();
            [KFem, KredFem] = comp.getElementLHS();

            comp.NewCreateIntegrators()
            KInt = comp.lhs;
            testCase.verifyEqual(KInt, KFem)
        end


        function testQuadrature(testCase, triangle)
            s.computerType    = 'FEM';
            s.testName         = triangle;
            s.variablesToStore = {'d_u'};
            test = PrecomputedVariableTest(s);
            comp = test.getComputation();
            elem = comp.getElement();
            quadFem = elem.quadrature;

            comp.NewCreateIntegrators()
            quadInt = comp.integrator.lhsint.getQuadrature;
            testCase.verifyEqual(quadInt, quadFem)
        end
    end


end