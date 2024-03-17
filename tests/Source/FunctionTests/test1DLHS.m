classdef test1DLHS < handle

    properties (Access = public)
        tol = 1e-12;
    end

    properties (Access = private)
        mesh
        K, M
        LHS, RHS
        funL2, funH1
    end

    methods (Access = protected, Static)
        
        function nu = createPoissonValue()
            nu = 0;
        end
    end
    
    methods (Access = public)
        
        function obj = test1DLHS()
            obj.init()
            obj.createMesh()
            obj.createLHS();
            obj.createRHS();
            obj.solveSystem();
        end

        function p = hasPassed(obj)
            x = load('test_1DLHS', 'x');
            p = isequal(x.x, obj.funH1.fValues);
        end
        
    end

    methods (Access = private)

        function init(obj)
        end

        function createMesh(obj)
            n = 50;
            xCoord =  0:1/(n-1):1;
            s.coord(:,1) = xCoord;
            s.coord(:,2) = zeros(n,1);
            s.connec(:,1) = 1:n-1;
            s.connec(:,2) = 2:n;
            s.kFace = -1;
            m = Mesh.create(s);
%             m.plot();
            obj.mesh = m;
        end

        function createLHS(obj)
            obj.createStifnessMatrix();
            obj.createMassMatrix();
            e = 0.1*obj.mesh.computeMeanCellSize();
            obj.LHS = e*obj.K + obj.M;
        end

        function createStifnessMatrix(obj)
            l.type  = 'StiffnessMatrix';
            l.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            l.test  = LagrangianFunction.create(obj.mesh,1,'P1');
            l.mesh  = obj.mesh;
            lhs = LHSintegrator.create(l);
            obj.K = lhs.compute();
        end
        
        function createMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.test  = LagrangianFunction.create(obj.mesh,1,'P1');
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh  = obj.mesh;
            lhs    = LHSintegrator.create(s);
            obj.M  = lhs.compute();
        end

        function createRHS(obj)
            xCoord = obj.mesh.coord(:,1)';
            xmean = (max(xCoord)-min(xCoord))/2;
            f = heaviside(xCoord-xmean);
            
            s.mesh = obj.mesh;
            s.fValues = f;
            s.order = 'P1';
            fL2 = LagrangianFunction(s);
%             fL2.plot;
            obj.RHS = obj.M*(fL2.fValues)';
            obj.funL2 = fL2;
        end

        function solveSystem(obj)
            s.type    =  'DIRECT';
            solver    = Solver.create(s);
            fH1Values = solver.solve(obj.LHS,obj.RHS);
    
            s.mesh = obj.mesh;
            s.fValues = fH1Values;
            s.order = 'P1';
            fH1 = LagrangianFunction(s);
%             fH1.plot()
            obj.funH1 = fH1;
%             obj.mesh.plot();
%             hold on
%             xCoord = obj.mesh.coord(:,1)';
%             plot(xCoord,obj.funL2.fValues,'-+');
%             plot(xCoord,obj.funH1.fValues,'-+');
        end
    end
    
end