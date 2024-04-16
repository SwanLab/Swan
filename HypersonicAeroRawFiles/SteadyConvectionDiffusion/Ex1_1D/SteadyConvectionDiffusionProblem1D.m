classdef SteadyConvectionDiffusionProblem1D < handle

% PER P1 FUNCIONA MILLOR QUE PER P2 !!!

% This program solves a one-dimensional convection-diffusion equation
%        a·u_x - nu·u_xx = f
% with Dirichlet boundary conditions using
% the finite element method with some stabilized formulations.

%  START INPUTS

% problem
% disp('This program solves a convection-diffusion equation')
% disp('         a·u_x - nu·u_xx = f')
% disp('with 0<x<1 and essential boundary conditions on both ends.')
% disp(' ')
% disp('One of the following problems can be solved at once:')
% disp('[1]:  Boundary conditions: u(0)= 0, u(1) = 1 ')
% disp('      Source term: f = 0')
% disp('[2]:  Boundary conditions: u(0)= 0, u(1) = 0 ')
% disp('      Source term: f = 1')
% disp('[3]:  Boundary conditions: u(0)= 0, u(1) = 1 ')
% disp('      Source term: f = sin(pi·x)')
% disp('[4]:  Boundary conditions: u(0)= 0, u(1) = 1 ')
% disp('      Source term: f = 20*exp(-5*(x-1/8))-10*exp(-5*(x-1/4))')
% disp('[5]:  Boundary conditions: u(0)= 0, u(1) = 1 ')
% disp('      Source term: f = 10*exp(-5*x)-4*exp(-x)')

% a Convection coefficient
% nu Diffusion coefficient
% 100 number of elements
% p 1 lineal; 2 parabolic elements

% stab
% disp ('The problem can be solved using one of the following methods: ');
% disp ('    [1] Galerkin')
% disp ('    [2] Streamline Upwind (SU)');
% disp ('    [3] Streamline Upwind Petrov-Galerkin (SUPG)');
% disp ('    [4] Galerkin Least-Squares (GLS)');
% disp ('    [5] Sub-Grid Scale (SGS)');

% END INPUTS

    properties (Access = private)
        problem
        numel
        p
        trial
        stab
    end

    properties (Access = private)
        source
        coord
        mesh
    end

    methods (Access = public)
        function obj = SteadyConvectionDiffusionProblem1D(cParams)
            obj.init(cParams);
            obj.computeCoordinateMatrix();
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
            obj.applyPostProcess(sol,a,nu);
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

        function computeCoordinateMatrix(obj)
            nEl           = obj.numel;
            xnode         = 0:1/nEl:1;
            obj.coord     = xnode;
            s.coord       = xnode';
            s.connec(:,1) = 1:length(xnode)-1;
            s.connec(:,2) = 2:length(xnode);
            obj.mesh      = Mesh.create(s);
        end

        function [K,f] = computeSystemLHSandRHS(obj,a,nu)
            h  = obj.coord(2)-obj.coord(1);
            Pe = a*h/(2*nu);
            if obj.stab==1
                s.trial = obj.trial;
                s.stab = obj.stab;
                s.mesh = obj.mesh;
                wf     = WeakFormSolver.create(s);
                [K,f]  = wf.compute(a,nu,obj.source);
            elseif obj.stab==2
                alfa = coth(Pe)-1/Pe;
                tau = alfa*h/(2*a); % Recommended
                if obj.p == 1
                    s.trial = obj.trial;
                    s.stab = obj.stab;
                    s.mesh = obj.mesh;
                    s.tau  = tau;
                    wf     = WeakFormSolver.create(s);
                    [K,f]  = wf.compute(a,nu,obj.source);
                else
                    beta = 2*((coth(Pe)-1/Pe)-(cosh(Pe))^2*(coth(2*Pe)-1/(2*Pe)))/(2-(cosh(Pe))^2);
                    tau_c = beta*h/(2*a); % Recommended
                    s.trial = obj.trial;
                    s.stab = obj.stab;
                    s.mesh = obj.mesh;
                    s.tau  = diag([tau_c,tau_c,tau]); % must be function
                    wf     = WeakFormSolver.create(s);
                    [K,f]  = wf.compute(a,nu,obj.source);
                end
            elseif obj.stab == 3
                alfa = coth(Pe)-1/Pe;
                tau = alfa*h/(2*a); % Recommended
                if obj.p == 1
                    [K,f] = system_SUPG_p1(tau,a,nu,obj.coord,obj.source);
                else
                    beta = ((2*Pe-1)*exp(3*Pe)+(-6*Pe+7)*exp(Pe)+(-6*Pe-7)*exp(-Pe)+(2*Pe+1)*exp(-3*Pe))...
                        /( (Pe+3)*exp(3*Pe)+(-7*Pe-3)*exp(Pe)+ (7*Pe-3)*exp(-Pe)-  (Pe+3)*exp(-3*Pe));
                    tau_c = beta*h/(2*a); % Recommended
                    [K,f] = system_SUPG_p2(tau, tau_c,a,nu,obj.coord,obj.source);
                end
            elseif obj.stab == 4
                alfa = coth(Pe)-1/Pe;
                tau = alfa*h/(2*a); % Recommended
                if obj.p == 1
                    % When using linear elements, the same solution is obtained with SUPG and GLS
                    [K,f] = system_SUPG_p1(tau,a,nu,obj.coord,obj.source);
                else
                    beta = ((2*Pe-1)*exp(3*Pe)+(-6*Pe+7)*exp(Pe)+(-6*Pe-7)*exp(-Pe)+(2*Pe+1)*exp(-3*Pe))...
                        /( (Pe+3)*exp(3*Pe)+(-7*Pe-3)*exp(Pe)+ (7*Pe-3)*exp(-Pe)-  (Pe+3)*exp(-3*Pe));
                    tau_c = beta*h/(2*a); % Recommended
                    [K,f] = system_GLS_p2(tau, tau_c,a,nu,obj.coord,obj.source);
                end
            elseif obj.stab == 5
                alfa = coth(Pe)-1/Pe;
                tau = alfa*h/(2*a); % Recommended
                if obj.p == 1
                    % When using linear elements, the same solution is obtained with SUPG and SGS
                    [K,f] = system_SUPG_p1(tau,a,nu,obj.coord,obj.source);
                else
                    beta = ((2*Pe-1)*exp(3*Pe)+(-6*Pe+7)*exp(Pe)+(-6*Pe-7)*exp(-Pe)+(2*Pe+1)*exp(-3*Pe))...
                        /( (Pe+3)*exp(3*Pe)+(-7*Pe-3)*exp(Pe)+ (7*Pe-3)*exp(-Pe)-  (Pe+3)*exp(-3*Pe));
                    tau_c = beta*h/(2*a); % Recommended
                    [K,f] = system_SGS_p2(tau, tau_c,a,nu,obj.coord,obj.source);
                end
            end
        end

        function sol = solveSystem(obj,K,f)
            s.nnodes    = length(obj.coord);
            s.p         = obj.p;
            s.K         = K;
            s.f         = f;
            KCf         = MonolithicFormComputer(s);
            [Ktot,ftot] = KCf.compute(obj.problem);
            sol         = Ktot\ftot;
            obj.trial.fValues = sol(1:end-2);
        end

        function applyPostProcess(obj,sol,a,nu)
            uh     = sol(1:end-2);
            x      = 0:0.01:1;
            exacto = u(x,a,nu,obj.problem);
            if obj.p == 1
                plot(obj.coord,uh,'r-',x,exacto,'k:','LineWidth',3,'MarkerSize',10)
            else
                [x0,y0]=obj.plotQuadraticElements();
                plot(x0,y0,'r-', x,exacto,'k:','LineWidth',3,'MarkerSize',10)
            end
            if obj.stab == 1
                l = legend('Galerkin solution','exact solution');
            elseif obj.stab==2
                l = legend('SU solution','exact solution');
            elseif obj.stab == 3
                l = legend('SUPG solution','exact solution');
            elseif obj.stab == 4
                l = legend('GLS solution','exact solution');
            else
                l = legend('SGS solution','exact solution');
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