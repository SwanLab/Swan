classdef DomainFunTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        cases = {'DDP','Grad','Partial','SymGrad','InternalEnergy','VolumetricEnergy','DeviatoricEnergy'}
    end

    methods (Test, TestTags = {'DomainFun', 'LevelSet'})
        function testLevelSet2D(testCase, cases)
            filename = ['testDomainFun',cases,'2D'];
            m        = testCase.obtain2DTestMesh();


            s.type   = 'Full';
            g        = GeometricalFunction(s);
            lsFun    = g.computeLevelSetFunction(m);
            s.fun  = lsFun;
            s.mesh = m;
            s.type = 'LevelSet';
            s.plotting = false;
            ls     = DesignVariable.create(s);



            xNew     = lsFun.fValues;
            load(filename,'xRef');
            err      = norm(xNew-xRef)/norm(xRef);
            tol      = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

        function testLevelSet3D(testCase, cases)
            filename = ['testDomainFun',cases,'3D'];
            m        = testCase.obtain3DTestMesh();
            s.type   = cases;
            g        = GeometricalFunction(s);
            lsFun    = g.computeLevelSetFunction(m);
            xNew     = lsFun.fValues;
            load(filename,'xRef');
            err      = norm(xNew-xRef)/norm(xRef);
            tol      = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end
    end

    methods (Static, Access = private)
        function m = obtain2DTestMesh()
            m = QuadMesh(2,1,100,50);
        end

        function m = obtain3DTestMesh()
            m = HexaMesh(2,1,1,20,10,10);
        end

    end

end