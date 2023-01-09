classdef RemeshingTestsPrueba < handle 
   
    methods (Access = public)

        function obj = RemeshingTestsPrueba()
           obj.testRemeshP1DiscontinousFunction();
        end

    end



    methods (Access = private)

        function testRemeshP1ContinousFunction(obj)
            mC = obj.createCoarseMesh();
            f  = obj.createP1ContinousFunction(mC);
            mF = mC.remesh();        
            fFine = f.refine(mC,mF);   
        end

        function testRemeshP1DiscontinousFunction(obj)
            mC = obj.createCoarseDiscontinousMesh();
            fFine  = obj.createP1DiscontinousFunction(mC);
            for i = 1:2
            mF = mC.remesh();         
            fFine = fFine.refine(mC,mF);        
            mC    = mF.createDiscontinuousMesh();
            fFine.plot(mC);
            end
        end        

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