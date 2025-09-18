classdef Connect < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        mesh
        test
        epsilon
        charFun
        boundaryConditions
        diffCoef
        massCoef
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)

        function obj = Connect()
            obj.createMesh();
            obj.test = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.epsilon = 1e-3;
            obj.createCharacteristicFunction();
            obj.createBoundaryConditions();
            obj.createDiffusionCoeficient();
            obj.createMassCoeficient();
            LHS = obj.createLHS();
            RHS = obj.createRHS();


            s.mesh = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            bcA = BCApplier(s);

            uFun = LagrangianFunction.create(obj.mesh,1,'P1');
            sS.type      = 'DIRECT';
            solver       = Solver.create(sS);
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solver     = solver;
            s.boundaryConditions = obj.boundaryConditions;
            s.BCApplier          = bcA;
            problemSolver    = ProblemSolver(s);   

            s.stiffness = LHS;
            s.forces    = RHS;  

            [u,~]       = problemSolver.solve(s);           
            %uSplit = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            uFun.setFValues(u);     
            plot(uFun)
            plot(uFun.*obj.charFun)
        end
        
        function createDiffusionCoeficient(obj)
            e   = obj.epsilon;
            x   = obj.charFun;
            k0  = e;
            k1  = 1-e;            
            obj.diffCoef = k0*(1-x)+k1*x;
        end

        function createMassCoeficient(obj)
            e   = obj.epsilon;
            x   = obj.charFun;
            m0  = 1-e;
            m1  = e;
            obj.massCoef = m0*(1-x)+m1*x;
        end

        function LHS = createLHS(obj)   
            k   = obj.diffCoef;
            m   = obj.massCoef;
            LHS = IntegrateLHS(@(u,v) DDP(Grad(v),k.*Grad(u)) + DP(v,m.*u),obj.test,obj.test,obj.mesh);
        end

        function RHS = createRHS(obj)
            x   = obj.charFun;
            m   = obj.massCoef;
            RHS = IntegrateRHS(@(v) DP(v,m.*(1-x)),obj.test,obj.mesh);
        end




    end
    
    methods (Access = private)
        

        function createMesh(obj)
            obj.mesh = UnitTriangleMesh(20,20);
        end

        function createCharacteristicFunction(obj)
            s.dim = obj.mesh.ndim;
            s.nHoles = [1 1];
            s.phases = [0 0]; 
            s.phiZero = 0.5;
            s.totalLengths = [1,1];
            s.type         = 'Holes';
            g              = GeometricalFunction(s);
            phiFun         = g.computeLevelSetFunction(obj.mesh);
            ls          = phiFun.fValues;
            lsInclusion = -ls;
            sU.backgroundMesh = obj.mesh;
            sU.boundaryMesh   = obj.mesh.createBoundaryMesh;
            uMesh             = UnfittedMesh(sU);
            uMesh.compute(lsInclusion);           
            cFun = CharacteristicFunction.create(uMesh);

            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = obj.test;
            f = Filter.create(s);

            obj.charFun = f.compute(cFun,2);
        end

        function createBoundaryConditions(obj)
            isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
            isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
            isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2))) < 1e-12);
            isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2))) < 1e-12);
            
           
            sDir{1}.domain    = @(coor) isTop(coor) | isLeft(coor) | isBottom(coor) | isRight(coor);
            sDir{1}.direction = 1;
            sDir{1}.value     = 0;
            sDir{1}.ndim      = 1;

            dirichletFun = DirichletCondition(obj.mesh, sDir{1});

            s.dirichletFun = dirichletFun;
            s.pointloadFun = [];
            s.periodicFun  = [];
            s.mesh = obj.mesh;
            obj.boundaryConditions = BoundaryConditions(s);            
        end
        
        
        
    end
    
end