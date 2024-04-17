classdef SteadyConvectionDiffusionProblem < handle

    properties (Access = private)
        problem
        numel
        p
        trial
        stab
    end

    properties (Access = private)
        source
        mesh
    end

    methods (Access = public)
        function obj = SteadyConvectionDiffusionProblem(cParams)
            obj.init(cParams);
            obj.computeMesh();
            switch obj.p
                case 1
                    obj.trial = LagrangianFunction.create(obj.mesh,1,'P1');
                case 2
                    obj.trial = LagrangianFunction.create(obj.mesh,1,'P2');
            end
            obj.computeSourceTerm();
        end

        function sol = compute(obj,a,nu)
            [K,f] = obj.computeSystemLHSandRHS(a,nu);
            sol   = obj.solveSystem(K,f);
            obj.plotSolution(sol);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.problem = cParams.problem;
            obj.numel   = cParams.numel;
            obj.p       = cParams.p;
            obj.stab    = cParams.stab;
        end

        function computeSourceTerm(obj)
            switch obj.problem
                case 1
                    s = AnalyticalFunction.create(@(x) zeros(size(x(1,:,:))),1,obj.mesh);
                case 2
                    s = AnalyticalFunction.create(@(x) ones(size(x(1,:,:))), 1, obj.mesh);
                case 3
                    s = AnalyticalFunction.create(@(x) sin(pi*x(1,:,:)), 1, obj.mesh);
                case 4
                    s = AnalyticalFunction.create(@(x) 20*exp(-5*(x(1,:,:)-1/8))-10*exp(-5*(x(1,:,:)-1/4)), 1, obj.mesh);
                case 5
                    s = AnalyticalFunction.create(@(x) 10*exp(-5*x(1,:,:))-4*exp(-x(1,:,:)), 1, obj.mesh);
            end
            obj.source = s;
        end

        function computeMesh(obj)
            nEl           = obj.numel;
            xnode         = 0:1/nEl:1;
            s.coord       = xnode';
            s.connec(:,1) = 1:length(xnode)-1;
            s.connec(:,2) = 2:length(xnode);
            obj.mesh      = Mesh.create(s);
        end

        function [K,f] = computeSystemLHSandRHS(obj,a,nu) % taus shown below are recommended, but people from lacan allow user-defined taus also
            h  = obj.mesh.computeMeanCellSize();
            Pe = a*h/(2*nu);
            if obj.stab==1
                s.trial = obj.trial;
                s.stab = obj.stab;
                s.mesh = obj.mesh;
                wf     = WeakFormSolver.create(s);
                [K,f]  = wf.compute(a,nu,obj.source);
            elseif obj.stab==2
                alfa = coth(Pe)-1/Pe;
                tau = alfa*h/(2*a);
                if obj.p == 1
                    s.trial = obj.trial;
                    s.stab = obj.stab;
                    s.mesh = obj.mesh;
                    s.tau  = tau;
                    wf     = WeakFormSolver.create(s);
                    [K,f]  = wf.compute(a,nu,obj.source);
                else
                    beta = 2*((coth(Pe)-1/Pe)-(cosh(Pe))^2*(coth(2*Pe)-1/(2*Pe)))/(2-(cosh(Pe))^2);
                    tau_c = beta*h/(2*a);
                    s.trial = obj.trial;
                    s.stab = obj.stab;
                    s.mesh = obj.mesh;
                    s.tau  = diag([tau_c,tau_c,tau]);
                    wf     = WeakFormSolver.create(s);
                    [K,f]  = wf.compute(a,nu,obj.source);
                end
            elseif obj.stab == 3
                alfa = coth(Pe)-1/Pe;
                tau  = alfa*h/(2*a);
                if obj.p == 1
                    s.trial = obj.trial;
                    s.stab = obj.stab;
                    s.mesh = obj.mesh;
                    s.tau  = tau;
                    wf     = WeakFormSolver.create(s);
                    [K,f]  = wf.compute(a,nu,obj.source);
                else
                    error('SUPG not available with quadratic elements')
                end
            end
        end

        function sol = solveSystem(obj,K,f)
            s.nnodes    = obj.mesh.nnodes;
            s.p         = obj.p;
            s.K         = K;
            s.f         = f;
            KCf         = MonolithicFormComputer(s);
            [Ktot,ftot] = KCf.compute(obj.problem);
            sol         = Ktot\ftot;
            obj.trial.fValues = sol(1:end-2);
        end

        function plotSolution(obj,sol)
            x  = obj.mesh.coord;
            uh = sol(1:end-2);
            if obj.p == 1
                plot(x,uh,'r-','LineWidth',1.5)
            else
                [x0,y0]=obj.plotQuadraticElements();
                plot(x0,y0,'r-','LineWidth',1.5)
            end
            if obj.stab == 1
                l = legend('Galerkin solution');
            elseif obj.stab==2
                l = legend('SU solution');
            elseif obj.stab == 3
                l = legend('SUPG solution');
            end
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
    end
end