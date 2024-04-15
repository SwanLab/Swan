function sol = computeSteadyConvectionDiffusion1D(problem,a,nu,numel,p,stab)

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

switch problem
    case 1
        source = SourceTermComputer.create(@(x) 0);
    case 2
        source = SourceTermComputer.create(@(x) 1);
    case 3
        source = SourceTermComputer.create(@(x) sin(pi*x));
    case 4
        source = SourceTermComputer.create(@(x) 20*exp(-5*(x-1/8))-10*exp(-5*(x-1/4)));
    case 5
        source = SourceTermComputer.create(@(x) 10*exp(-5*x)-4*exp(-x));
end

if p == 1 
    h = 1/numel;
    numnp = numel + 1;
elseif p == 2
    h = 1/(2*numel);
    numnp = 2*numel + 1;
else
    error('Unavaible elements')
end
xnode = 0:h:1; 
  
Pe = a*h/(2*nu);
disp(' ')
disp(strcat('Peclet number: ',num2str(Pe)));
 
% Stabilizatin parameters and system obtained after discretization
if stab==1
    disp('Galerkin formulation')
    if p==1
        [K,f] = system_p1(a,nu,xnode,source); 
    else
        [K,f] = system_p2(a,nu,xnode,source);       
    end
elseif stab==2 
    disp('SU formulation')
    alfa = coth(Pe)-1/Pe; 
    tau = alfa*h/(2*a); % Recommended
    if isempty(tau)
        tau = alfa*h/(2*a);
    end
    if p == 1
        [K,f] = system_SU_p1(tau,a,nu,xnode,source); 
    else
        beta = 2*((coth(Pe)-1/Pe)-(cosh(Pe))^2*(coth(2*Pe)-1/(2*Pe)))/(2-(cosh(Pe))^2);
        tau_c = beta*h/(2*a); % Recommended
        if isempty(tau_c)
            tau_c = beta*h/(2*a);
        end
        [K,f] = system_SU_p2(tau, tau_c,a,nu,xnode,source);
    end
elseif stab == 3
    disp('SUPG formulation')
    alfa = coth(Pe)-1/Pe; 
    tau = alfa*h/(2*a); % Recommended
    if isempty(tau)
        tau = alfa*h/(2*a);
    end
    if p == 1 
        [K,f] = system_SUPG_p1(tau,a,nu,xnode,source); 
    else
        beta = ((2*Pe-1)*exp(3*Pe)+(-6*Pe+7)*exp(Pe)+(-6*Pe-7)*exp(-Pe)+(2*Pe+1)*exp(-3*Pe))...
               /( (Pe+3)*exp(3*Pe)+(-7*Pe-3)*exp(Pe)+ (7*Pe-3)*exp(-Pe)-  (Pe+3)*exp(-3*Pe));
        tau_c = beta*h/(2*a); % Recommended
        if isempty(tau_c)
            tau_c = beta*h/(2*a);
        end
        [K,f] = system_SUPG_p2(tau, tau_c,a,nu,xnode,source);
    end
elseif stab == 4
    disp ('GLS formulation');
    alfa = coth(Pe)-1/Pe; 
    tau = alfa*h/(2*a); % Recommended
    if isempty(tau)
        tau = alfa*h/(2*a);
    end
    if p == 1 
        % When using linear elements, the same solution is obtained with SUPG and GLS 
        [K,f] = system_SUPG_p1(tau,a,nu,xnode,source); 
    else
        beta = ((2*Pe-1)*exp(3*Pe)+(-6*Pe+7)*exp(Pe)+(-6*Pe-7)*exp(-Pe)+(2*Pe+1)*exp(-3*Pe))...
               /( (Pe+3)*exp(3*Pe)+(-7*Pe-3)*exp(Pe)+ (7*Pe-3)*exp(-Pe)-  (Pe+3)*exp(-3*Pe));
        tau_c = beta*h/(2*a); % Recommended
        if isempty(tau_c)
            tau_c = beta*h/(2*a);
        end
        [K,f] = system_GLS_p2(tau, tau_c,a,nu,xnode,source);
    end
elseif stab == 5
    disp ('Formulación SGS');
    alfa = coth(Pe)-1/Pe; 
    tau = alfa*h/(2*a); % Recommended
    if isempty(tau)
        tau = alfa*h/(2*a);
    end
    if p == 1
        % When using linear elements, the same solution is obtained with SUPG and SGS 
        [K,f] = system_SUPG_p1(tau,a,nu,xnode,source); 
    else
        beta = ((2*Pe-1)*exp(3*Pe)+(-6*Pe+7)*exp(Pe)+(-6*Pe-7)*exp(-Pe)+(2*Pe+1)*exp(-3*Pe))...
               /( (Pe+3)*exp(3*Pe)+(-7*Pe-3)*exp(Pe)+ (7*Pe-3)*exp(-Pe)-  (Pe+3)*exp(-3*Pe));
        tau_c = beta*h/(2*a); % Recommended
        if isempty(tau_c)
            tau_c = beta*h/(2*a);
        end
        [K,f] = system_SGS_p2(tau, tau_c,a,nu,xnode,source);
    end
else
    error('unavailable formulation');
end

s.nnodes    = numnp;
s.K         = K;
s.f         = f;
KCf         = MonolithicFormComputer(s);
[Ktot,ftot] = KCf.compute(problem);

% Solution 
sol = Ktot\ftot;
  
% Postproces
u_h = sol(1:numnp);
x = 0:0.01:1;
exacto = u(x,a,nu,problem);
if p == 1
    plot(xnode,u_h,'r-o',x,exacto,'k:','LineWidth',3,'MarkerSize',10)
else
    [x0,y0]=plot_p2(sol,xnode);
    plot(xnode,u_h,'ro', x,exacto,'k:', x0,y0,'r-','LineWidth',3,'MarkerSize',10)
end
if stab == 1
    l = legend('Galerkin solution','exact solution');
elseif stab==2
    l = legend('SU solution','exact solution');
elseif stab == 3
    l = legend('SUPG solution','exact solution');
elseif stab == 4
    l = legend('GLS solution','exact solution');
else   
    l = legend('SGS solution','exact solution');
end
set(l, 'FontSize',14);
set(gca, 'FontSize',14);
 
end