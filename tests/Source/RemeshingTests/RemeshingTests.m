classdef RemeshingTests < handle & matlab.unittest.TestCase


    methods (Test, TestTags = {'Remesh'})

        function testRemeshMesh(obj)
            mC = obj.createCoarseMesh();
            mF = mC.remesh(1);
            s = load('test_RemeshMesh');
            err(1) = norm(mF.coord(:)  - s.meshFine.coord(:));
            err(2) = norm(mF.connec(:) - s.meshFine.connec(:));
            tol = 1e-6;
            obj.verifyLessThanOrEqual(norm(err), tol)
        end

        function testRemeshP1ContinousFunction(obj)
            m = obj.createCoarseMesh();
            f = obj.createP1ContinousFunction(m);
            for i = 1:2
                mF = m.remesh(1);
                f = f.refine(mF);
                m    = mF;
            end
            s = load('test_RemeshP1Function');
            err = norm(s.fValues(:) - f.fValues(:));
            tol = 1e-6;
            obj.verifyLessThanOrEqual(err, tol)
        end

        function testRemeshP1DiscontinousFunction(obj)
            m = obj.createCoarseDiscontinousMesh();
            f = obj.createP1DiscontinousFunction(m);
            for i = 1:2
                mF = m.remesh(1); % mF is CONTINUOUS, m is DISCONTINUOUS
                f = f.refine(m); %fNew = fOld.refine(mF); -> fOld has m, fNew has mF
                m = mF.createDiscontinuousMesh(); % care: discont/cont -> interpolation through continuous
                % ideally: work only with discontinuous meshes
            end
            s = load('test_RemeshP1DiscFunction');
            err = norm(s.fValues(:) - f.fValues(:));
            tol = 1e-6;
            obj.verifyLessThanOrEqual(err, tol)
        end

    end

    methods (Access = private)

        function m = createCoarseMesh(obj)
            s.coord = [1 0; 0 1; 0 0; 1 1];
            s.connec = [3 1 2; 1 4 2];
            m = Mesh.create(s);
        end

        function mD = createCoarseDiscontinousMesh(obj)
            m   = obj.createCoarseMesh();
            mD = m.createDiscontinuousMesh();
        end

        function fC = createP1ContinousFunction(obj,m)
            f         = obj.createFunctionToRemesh();
            s.mesh    = m;
            s.fValues = f(m.coord);
            s.order   = 'P1';
            fC        = LagrangianFunction(s);
        end

        function fC = createP1DiscontinousFunction(obj,m)
            f         = obj.createFunctionToRemesh();
            s.fValues = f(m.computeBaricenter()');
            s.mesh    = m;
            s.order   = 'P0';
            f0 = LagrangianFunction(s);
            fC = f0.project('P1D');
        end

        function f = createFunctionToRemesh(obj)
            f = @(x) x(:,1).^2+x(:,2).^2;
        end

    end

end