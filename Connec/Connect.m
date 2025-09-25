classdef Connect < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        mesh
        diffCoef
        massCoef
        test
        boundaryConditions
        problemSolver
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)

        function obj = Connect(cParams)
            obj.init(cParams);
            obj.test = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.createBoundaryConditions();
            obj.createSolver();
        end
        
        function [uFun] = solve(obj,x)
            k   = obj.diffCoef(x);
            m   = obj.massCoef(x);
            LHS = IntegrateLHS(@(u,v) DDP(Grad(v),k.*Grad(u)) + DP(v,m.*u),obj.test,obj.test,obj.mesh);
            RHS = IntegrateRHS(@(v) DP(v,m.*(x)),obj.test,obj.mesh);
            s.stiffness = LHS;
            s.forces    = RHS;  
            [u,~]       = obj.problemSolver.solve(s); 
            uFun = LagrangianFunction.create(obj.mesh,1,'P1');
            uFun.setFValues(u);  
        end

    end
    
    methods (Access = private)
        

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.diffCoef = cParams.diffCoef;
            obj.massCoef = cParams.massCoef;
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

        function createSolver(obj)
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;           
            bcA = BCApplier(s);
            sS.type      = 'DIRECT';
            solver       = Solver.create(sS);
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solver     = solver;
            s.boundaryConditions = obj.boundaryConditions;
            s.BCApplier          = bcA;
            obj.problemSolver    = ProblemSolver(s);   
       end             
        
        
    end
    
end