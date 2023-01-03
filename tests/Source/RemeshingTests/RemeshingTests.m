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
            mC = obj.createCoarseMesh();
            f  = obj.createP1ContinousFunction(mC);
            mF = mC.remesh();        
            fFine = f.refine(mC,mF);   
            s = load('test_RemeshP1Function');
            err = norm(s.fValues(:) - fFine.fValues(:));
            tol = 1e-6;
            obj.verifyLessThanOrEqual(err, tol)
        end

        function testRemeshP1DiscontinousFunction(obj)
            mC = obj.createCoarseDiscontinousMesh();
            f  = obj.createP1DiscontinousFunction(mC);
            mF = mC.remesh();        
            fFine = f.refine(mC,mF);   
            s = load('test_RemeshP1DiscFunction');
            err = norm(s.fValues(:) - fFine.fValues(:));
            tol = 1e-6;
            obj.verifyLessThanOrEqual(err, tol)
        end        

    end

    methods (Access = private)

        function m = createCoarseMesh(obj)
            s.coord = [1 0; 0 1; 0 0; 1 1];
            s.connec = [3 1 2; 1 4 2];
            m = Mesh(s);
        end     

        function mD = createCoarseDiscontinousMesh(obj)
            m   = obj.createCoarseMesh;
            mD = m.createDiscontinuousMesh();
        end                

        function fC = createP1ContinousFunction(obj,m)               
            f         = obj.createFunctionToRemesh();
            s.type    = m.type;
            s.connec  = m.connec;
            s.fValues = f(m.coord);    
            fC        = P1Function(s);
        end    

        function fC = createP1DiscontinousFunction(obj,m)
            f         = obj.createFunctionToRemesh();
            s.fValues = f(m.computeBaricenter()');
            s.connec  = m.connec;
            s.type    = m.type;
            f0 = P0Function(s);
            fC = f0.computeP1DiscontinuousFunction(m);
        end        

        function f = createFunctionToRemesh(obj)
            f = @(x) x(:,1).*x(:,2);
        end        

    end

end