classdef ConvectionDiffusionTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        tests  = {'lowPecletLinearElem','highPecletQuadElem'};
        p1test = {'lowPecletLinearElem'};
    end

    methods (Test, TestTags = {'ConvectionDiffusion'})
        function testWithoutSourceTermGalerkin(testCase, tests)
            run(tests)
            nEl            = 100;
            xnode          = 0:1/nEl:1;
            sM.coord       = xnode';
            sM.connec(:,1) = 1:length(xnode)-1;
            sM.connec(:,2) = 2:length(xnode);
            m              = Mesh.create(sM);
            switch p
                case 1
                    s.trial = LagrangianFunction.create(m,1,'P1');
                case 2
                    s.trial = LagrangianFunction.create(m,1,'P2');
            end
            s.mesh    = m;
            s.sHandle = @(x) zeros(size(x(1,:,:)));
            s.dirValues.left  = 0;
            s.dirValues.right = 1;
            s.numel   = 100;
            s.stab    = 'Galerkin';
            prob      = SteadyConvectionDiffusionProblem(s);
            prob.compute(a,nu);
            sol = prob.trial.fValues;
            load([tests,'Steady1DGalerkinNoSource.mat'],'x');
            x = x(1:end-2);
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

        function testWithSourceTermGalerkin(testCase,tests)
            run(tests)
            nEl            = 100;
            xnode          = 0:1/nEl:1;
            sM.coord       = xnode';
            sM.connec(:,1) = 1:length(xnode)-1;
            sM.connec(:,2) = 2:length(xnode);
            m              = Mesh.create(sM);
            switch p
                case 1
                    s.trial = LagrangianFunction.create(m,1,'P1');
                case 2
                    s.trial = LagrangianFunction.create(m,1,'P2');
            end
            s.mesh    = m;
            s.sHandle = @(x) sin(pi*x(1,:,:));
            s.dirValues.left  = 0;
            s.dirValues.right = 1;
            s.numel   = 100;
            s.stab    = 'Galerkin';
            prob      = SteadyConvectionDiffusionProblem(s);
            prob.compute(a,nu);
            sol = prob.trial.fValues;
            load([tests,'Steady1DGalerkinSource.mat'],'x');
            x = x(1:end-2);
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

        function testWithoutSourceTermSU(testCase,tests)
            run(tests)
            nEl            = 100;
            xnode          = 0:1/nEl:1;
            sM.coord       = xnode';
            sM.connec(:,1) = 1:length(xnode)-1;
            sM.connec(:,2) = 2:length(xnode);
            m              = Mesh.create(sM);
            switch p
                case 1
                    s.trial = LagrangianFunction.create(m,1,'P1');
                case 2
                    s.trial = LagrangianFunction.create(m,1,'P2');
            end
            s.mesh    = m;
            s.sHandle = @(x) zeros(size(x(1,:,:)));
            s.dirValues.left  = 0;
            s.dirValues.right = 1;
            s.numel   = 100;
            s.stab    = 'Upwind';
            prob      = SteadyConvectionDiffusionProblem(s);
            prob.compute(a,nu);
            sol = prob.trial.fValues;
            load([tests,'Steady1DSUNoSource.mat'],'x');
            x = x(1:end-2);
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

        function testWithSourceTermSU(testCase,tests)
            run(tests)
            nEl            = 100;
            xnode          = 0:1/nEl:1;
            sM.coord       = xnode';
            sM.connec(:,1) = 1:length(xnode)-1;
            sM.connec(:,2) = 2:length(xnode);
            m              = Mesh.create(sM);
            switch p
                case 1
                    s.trial = LagrangianFunction.create(m,1,'P1');
                case 2
                    s.trial = LagrangianFunction.create(m,1,'P2');
            end
            s.mesh    = m;
            s.sHandle = @(x) sin(pi*x(1,:,:));
            s.dirValues.left  = 0;
            s.dirValues.right = 1;
            s.numel   = 100;
            s.stab    = 'Upwind';
            prob      = SteadyConvectionDiffusionProblem(s);
            prob.compute(a,nu);
            sol = prob.trial.fValues;
            load([tests,'Steady1DSUSource.mat'],'x');
            x = x(1:end-2);
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

        function testWithoutSourceTermSUPG(testCase,p1test)
            run(p1test)
            nEl            = 100;
            xnode          = 0:1/nEl:1;
            sM.coord       = xnode';
            sM.connec(:,1) = 1:length(xnode)-1;
            sM.connec(:,2) = 2:length(xnode);
            m              = Mesh.create(sM);
            switch p
                case 1
                    s.trial = LagrangianFunction.create(m,1,'P1');
                case 2
                    s.trial = LagrangianFunction.create(m,1,'P2');
            end
            s.mesh    = m;
            s.sHandle = @(x) zeros(size(x(1,:,:)));
            s.dirValues.left  = 0;
            s.dirValues.right = 1;
            s.numel   = 100;
            s.stab    = 'SUPG';
            prob      = SteadyConvectionDiffusionProblem(s);
            prob.compute(a,nu);
            sol = prob.trial.fValues;
            load([p1test,'Steady1DSUPGNoSource.mat'],'x');
            x = x(1:end-2);
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

        function testWithSourceTermSUPG(testCase,p1test)
            run(p1test)
            nEl            = 100;
            xnode          = 0:1/nEl:1;
            sM.coord       = xnode';
            sM.connec(:,1) = 1:length(xnode)-1;
            sM.connec(:,2) = 2:length(xnode);
            m              = Mesh.create(sM);
            switch p
                case 1
                    s.trial = LagrangianFunction.create(m,1,'P1');
                case 2
                    s.trial = LagrangianFunction.create(m,1,'P2');
            end
            s.mesh    = m;
            s.sHandle = @(x) sin(pi*x(1,:,:));
            s.dirValues.left  = 0;
            s.dirValues.right = 1;
            s.numel   = 100;
            s.stab    = 'SUPG';
            prob      = SteadyConvectionDiffusionProblem(s);
            prob.compute(a,nu);
            sol = prob.trial.fValues;
            load([p1test,'Steady1DSUPGSource.mat'],'x');
            x = x(1:end-2);
            err = norm(sol-x)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all;
        end

    end
end