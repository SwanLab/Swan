classdef SteadyConvectionDiffusionProblem < handle

    properties (Access = private)
        dirValues
        mesh
        trial
        stab
    end

    properties (Access = private)
        source
    end

    methods (Access = public)
        function obj = SteadyConvectionDiffusionProblem(cParams)
            obj.init(cParams);
            obj.computeSourceTerm(cParams);
        end

        function sol = compute(obj,a,nu)
            [K,f] = obj.computeSystemLHSandRHS(a,nu);
            sol   = obj.solveSystem(K,f);
            obj.plotSolution(sol);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.dirValues = cParams.dirValues;
            obj.mesh      = cParams.mesh;
            obj.trial     = cParams.trial;
            obj.stab      = cParams.stab;
        end

        function computeSourceTerm(obj,cParams)
            sH         = cParams.sHandle;
            s          = AnalyticalFunction.create(sH,1,obj.mesh);
            obj.source = s;
        end

        function [K,f] = computeSystemLHSandRHS(obj,a,nu)
            tau     = obj.computeRecommendedStabilizationParameter(a,nu);
            s.trial = obj.trial;
            s.stab  = obj.stab; % ['LHSintegratior',obj.stab]
            s.mesh  = obj.mesh;
            s.tau   = tau;
            wf      = WeakFormSolver.create(s); % Inside the a-vector would change in 2D
            [K,f]   = wf.compute(a,nu,obj.source);
        end

        function tau = computeRecommendedStabilizationParameter(obj,a,nu) % Independent class to choose tau between 1D/2D
            h  = obj.mesh.computeMeanCellSize(); % TAU AS LAGRANG FUN FIELD + include in LHSIntegrators
            Pe = a*h/(2*nu);
            switch obj.stab
                case 'Galerkin'
                    tau = [];
                case 'Upwind'
                    alfa = coth(Pe)-1/Pe;
                    tau  = alfa*h/(2*a);
                    if obj.trial.order == "P2"
                        beta = 2*((coth(Pe)-1/Pe)-(cosh(Pe))^2*(coth(2*Pe)-1/(2*Pe)))/(2-(cosh(Pe))^2);
                        tau_c = beta*h/(2*a);
                        tau  = diag([tau_c,tau_c,tau]);
                    end
                case 'SUPG'
                    alfa = coth(Pe)-1/Pe;
                    tau  = alfa*h/(2*a);
                    if obj.trial.order == "P2"
                        error('SUPG not available with quadratic elements')
                    end
            end
        end

        function sol = solveSystem(obj,K,f) % Fer aquí
            s.nnodes    = obj.mesh.nnodes;
            s.order     = obj.trial.order;
            s.K         = K;
            s.f         = f;
            s.bc        = obj.createBoundaryConditions1D();
            KCf         = MonolithicFormComputer(s);
            [Ktot,ftot] = KCf.compute();
            sol         = Ktot\ftot;
            obj.trial.fValues = sol(1:end-2);
        end

        function plotSolution(obj,sol) % Should differentiate between 1D and 2D
            x  = obj.mesh.coord;
            uh = sol(1:end-2);
            if obj.trial.order == "P1"
                plot(x,uh,'-','LineWidth',1.5)
            else
                [x0,y0]=obj.plotQuadraticElements();
                plot(x0,y0,'-','LineWidth',1.5)
            end
            l = legend([obj.stab,' solution']);
            set(l, 'FontSize',14);
            set(gca, 'FontSize',14);
        end

        function [x0,y0] = plotQuadraticElements(obj)
            xg     = -1:0.1:1;
            xElg   = obj.mesh.computeXgauss(xg);
            yElg   = obj.trial.evaluate(xg);
            [x0,v] = unique(xElg);
            y0     = yElg(v);
        end

        function bc = createBoundaryConditions1D(obj)
            xMin       = min(obj.mesh.coord(:,1));
            xMax       = max(obj.mesh.coord(:,1));
            isDirLeft  = @(coor)  abs(coor(:,1))==xMin;
            isDirRight = @(coor)  abs(coor(:,1))==xMax;

            sDir{1}.domain    = @(coor) isDirLeft(coor);
            sDir{1}.direction = 1;
            sDir{1}.value     = obj.dirValues.left;

            sDir{2}.domain    = @(coor) isDirRight(coor);
            sDir{2}.direction = 1;
            sDir{2}.value     = obj.dirValues.right;

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;
            s.pointloadFun = [];
            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);
        end
    end
end