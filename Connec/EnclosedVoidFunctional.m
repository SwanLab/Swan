classdef EnclosedVoidFunctional < handle

    properties (Access = private)
        value0
    end

    properties (Access = private)
        mesh
        filter
        boundaryConditions
        boundaryConditionsAdjont
        k
        m
        dm
        dk
        test
        problemSolver        
        problemSolverAdjoint
        target0
        target
        valueOld
    end

    methods (Access = public)
        function obj = EnclosedVoidFunctional(cParams)
            obj.init(cParams);
            obj.target = 0.01;
            obj.test = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.createBoundaryConditionsAdjoint();
            obj.createSolverAdjoint();
            obj.value0 = 1;
            obj.valueOld = -1;
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            obj.createBoundaryConditions(x);
            obj.createSolver();            
            xD  = x.obtainDomainFunction();           
            xR = obj.filterFields(xD);
            xR = xR{1};
            phi = obj.solveState(xR);%xR
            lam = obj.solveAdjoint(xR);
         %   rhoV = (phi - xD{1});
            rhoV = phi.*(1-xR);
            %rhoV = (1-phi).*massCoef(xR);
          %  plot(x.fun)
          %  plot(phi)
          %  plot(rhoV) 
          %  plot(lam);
            Dom   = Integrator.compute(ConstantFunction.create(1,obj.mesh),obj.mesh,2);
            J     = Integrator.compute(rhoV,obj.mesh,2)/Dom;

            
            %obj.updateEpsilonForNextIteration(J);
            %dJ{1} = -phi +  DP(obj.dm(xR).*(phi-xR)-obj.m(xR),lam) + 1.*obj.dk(xR).*DP(Grad(lam),Grad(phi));
            dJ{1} = -phi;% +  DP(obj.dm(xR).*(phi-xR)-obj.m(xR),lam) + 1.*obj.dk(xR).*DP(Grad(lam),Grad(phi));
            %dJ{1} =  -phi.*(1-xR);
            %dJ{1} =  -phi;%.*(1-xR);
            dJ{1} = dJ{1}./Dom*1;

            dJ = obj.filterFields(dJ);

            J  = obj.computeFunction(J);
            dJ = obj.computeGradient(dJ);
            %plot(dJ{1})
            %[J,dJ] = obj.computeComplianceFunctionAndGradient(x);
            obj.updateTarget0ForNextIteration(J);
        end


       function J = computeFunction(obj,P)
            pTar = obj.target0;
            J    = P/(pTar/obj.value0) - 1; % P-pTar/obj.value0 if pTar is close to zero!!
        end

        function dJ = computeGradient(obj,dP)
            pTar = obj.target0;
            dJ   = dP;
            dJ{1}.setFValues(dP{1}.fValues/(pTar/obj.value0));
        end

        function updateTarget0ForNextIteration(obj,J) 
            %if abs(J)<=1e-2
            if J-obj.valueOld<0 || abs(J) <= 1e-2              
               obj.target0 = max(obj.target0*(1-0.05),obj.target);
            end
            obj.valueOld = J;
        end        

    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.mesh     = cParams.mesh;
            obj.filter   = cParams.filter;
            obj.k = cParams.diffCoef;
            obj.m = cParams.massCoef;
            obj.dm       = cParams.dm;
            obj.dk       = cParams.dk;
            obj.target0  = cParams.target0;
            if isfield(cParams,'value0')
                obj.value0 = cParams.value0;
            end
        end

        function createBoundaryConditions(obj,x)
            [~, l2g]  = obj.mesh.createSingleBoundaryMesh();

            l2g = unique(l2g);
            fValues = x.fun.fValues(l2g);
           % xB = x.fun.restrictBaseToBoundary(bMesh);

            bF = x.fun.restrictToBoundary();
            isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
            isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
            isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2))) < 1e-12);
            isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2))) < 1e-12);
            sDir{1}.domain    = @(coor) isTop(coor) | isLeft(coor) | isBottom(coor) | isRight(coor);
            sDir{1}.direction = 1;
            sDir{1}.value     = fValues; %bF.fValues
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


        function xR = filterFields(obj,x)
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.filter.compute(x{i},2);
            end
        end        

        function [uFun] = solveState(obj,x)
            k   = obj.k(x);
            m   = obj.m(x);
            LHS = IntegrateLHS(@(u,v) DDP(Grad(v),k.*Grad(u)) + DP(v,m.*u),obj.test,obj.test,obj.mesh,'Domain');
            RHS = IntegrateRHS(@(v) DP(v,m.*x),obj.test,obj.mesh,'Domain');


            %LHS2 = IntegrateLHS(@(u,v) DP(v,u),obj.test,obj.test,obj.mesh,'Boundary');
            %RHS2 = IntegrateRHS(@(v) DP(v,x),obj.test,obj.mesh,'Boundary');


            s.stiffness = LHS;
            s.forces    = RHS;  
            [u,~]       = obj.problemSolver.solve(s); 
            uFun = LagrangianFunction.create(obj.mesh,1,'P1');
            uFun.setFValues(u);  
        end        

        function lam = solveAdjoint(obj,x)
            k   = obj.k(x);
            m   = obj.m(x);
            LHS = IntegrateLHS(@(u,v) DDP(Grad(v),k.*Grad(u)) + DP(v,m.*u),obj.test,obj.test,obj.mesh,'Domain');
            RHS = IntegrateRHS(@(v) -DP(v,m.*(1-x)),obj.test,obj.mesh,'Domain');
            s.stiffness = LHS;
            s.forces    = RHS;  
            [u,~]       = obj.problemSolverAdjoint.solve(s); 
            lam = LagrangianFunction.create(obj.mesh,1,'P1');
            lam.setFValues(u);
           % plot(1-x)
           % plot(-lam)
           % plot(-lam.*(1-x))
        end

        function createBoundaryConditionsAdjoint(obj)
            %[~, l2g]  = obj.mesh.createSingleBoundaryMesh();

            %l2g = unique(l2g);
            %fValues = x.fun.fValues(l2g);
           % xB = x.fun.restrictBaseToBoundary(bMesh);

            %bF = x.fun.restrictToBoundary();
            isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
            isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
            isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2))) < 1e-12);
            isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2))) < 1e-12);
            sDir{1}.domain    = @(coor) isTop(coor) | isLeft(coor) | isBottom(coor) | isRight(coor);
            sDir{1}.direction = 1;
            sDir{1}.value     = 0; %bF.fValues
            sDir{1}.ndim      = 1;
            dirichletFun = DirichletCondition(obj.mesh, sDir{1});
            s.dirichletFun = dirichletFun;
            s.pointloadFun = [];
            s.periodicFun  = [];
            s.mesh = obj.mesh;
            obj.boundaryConditionsAdjont = BoundaryConditions(s);            
        end        

        function createSolverAdjoint(obj)
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.boundaryConditionsAdjont;           
            bcA = BCApplier(s);
            sS.type      = 'DIRECT';
            solver       = Solver.create(sS);
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solver     = solver;
            s.boundaryConditions = obj.boundaryConditionsAdjont;
            s.BCApplier          = bcA;
            obj.problemSolverAdjoint    = ProblemSolver(s);   
        end          

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'EnclosedVolume';
        end
    end
end