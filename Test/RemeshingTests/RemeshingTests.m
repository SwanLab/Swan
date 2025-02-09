classdef RemeshingTests < handle & matlab.unittest.TestCase


    methods (Test, TestTags = {'Remesh'})

        function testRemeshMesh(obj)
            mC = obj.createCoarseMesh();
            mF = mC.remesh();
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
                m = m.remesh();
                f = f.refine(m);
            end
            s = load('test_RemeshP1Function');
            err = norm(s.fValues(:) - f.fValues(:));
            tol = 1e-6;
            obj.verifyLessThanOrEqual(err, tol)
        end

        function testRemeshP1DiscontinousFunction(obj)
            m = obj.createCoarseMesh();
            f = obj.createP1DiscontinousFunction(m);
            for i = 1:2
                m = m.remesh(); 
                f = f.refine(m); 
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

        function fC = createP1ContinousFunction(obj,m)
            f  = obj.createFunctionToRemesh(m);
            fC = f.project('P1');
        end

        function fC = createP1DiscontinousFunction(obj,m)
            f  = obj.createFunctionToRemesh(m);
            fC = f.project('P1D');
        end

        function f = createFunctionToRemesh(obj,mesh)
            s.fHandle = @(x) obj.parabola(x);
            s.mesh   = mesh;
            s.ndimf  = 1;
            f = AnalyticalFunction(s);
        end

        function f = parabola(obj,x)
            x1 = x(1,:,:);
            x2 = x(2,:,:);
            f  = x1.^2+x2.^2;
        end        


    end

end