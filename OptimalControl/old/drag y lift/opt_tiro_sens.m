function opt_tiro_sens
%-------------------------------------------------------------
%  ÓPTIMO DE TIRO PARABÓLICO CON DRAG & LIFT
%  Derivadas analíticas (sensibilidades variacionales)
%-------------------------------------------------------------
clc;  close all;

%% ------------------  PARÁMETROS FÍSICOS  -------------------
pars.g       = 9.81;                        % gravedad [m/s^2]
pars.rho     = 1.225;                       % densidad aire [kg/m^3]
pars.Sw      = 0.6;                         % superficie alar [m^2]
pars.Cd0     = 0.03;
pars.k       = 0.05;
pars.m       = 5;                           % masa [kg]
pars.Cl0     = 0;
pars.Clalpha = 0.1*360/(2*pi);              % [1/rad]  (0.1 por grado)

%% ------------------  DISCRETIZACIÓN -------------------------
N       = 200;                              % nº nodos de alfa(t)
pars.N  = N;
pars.s_nodes = linspace(0,1,N);             % malla normalizada 0…1
nu      = N + 1;                            % nº variables de decisión

%% ------------------  CONDICIONES INICIALES -----------------
x1_0   = 0;  x2_0 = 0;  v0 = 15;
gamma0 = deg2rad(30);                       % ángulo de disparo
y0     = [x1_0; x2_0; v0; gamma0];

%% ------------------  VARIABLES DE DECISIÓN -----------------
tf0      = 10;                              % conjetura tf
alpha0   = deg2rad(4);                      % conjetura α(t)
u0       = [ tf0,  ones(1,N)*alpha0 ];

lb       = [ 0.1,  ones(1,N)*deg2rad(-10) ];
ub       = [ 15 ,  ones(1,N)*deg2rad( 10) ];
u0       = max(lb, min(ub,u0));             % clip → dentro de límites

%% ------------------  OPCIONES FMINCON ----------------------
options = optimoptions('fmincon', ...
    'Algorithm',                 'sqp', ...
    'SpecifyObjectiveGradient',  true, ...
    'SpecifyConstraintGradient', true, ...
    'UseParallel',               true, ...
    'FiniteDifferenceType',      'central', ...   % sólo por si falláramos
    'Display',                   'iter');

%% ------------------  OPTIMIZACIÓN --------------------------
costfun  = @(u) cost_with_grad   (u,y0,pars);
confun   = @(u) con_with_grad    (u,y0,pars);

[u_opt,~,exitflag,output] = ...
    fmincon(costfun,u0,[],[],[],[],lb,ub,confun,options);

%% ------------------  RESULTADOS ----------------------------
fprintf('\n·  Alcance máximo    = %.3f  m\n', -costfun(u_opt));
fprintf('·  Tiempo de vuelo   = %.3f  s\n',  u_opt(1));
fprintf('·  Iteraciones       = %d\n',      output.iterations);
fprintf('·  Mensaje solver    : %s\n\n',    output.message);

% Trayectoria óptima para visualizar
[Y,~] = simulate_state_sens(u_opt,y0,pars);
figure;  plot(Y(:,1),Y(:,2),'LineWidth',1.6);
xlabel('Distancia horizontal  x  [m]');
ylabel('Altura  y  [m]');  grid on;  title('Trayectoria óptima');
end
%% =====================================================================
%                       COSTE  +  GRADIENTE
% ======================================================================
function [J, gradJ] = cost_with_grad(u,y0,pars)
[Y,S] = simulate_state_sens(u,y0,pars);

J         = -Y(end,1);            % max alcance → min ( -x(tf) )
gradJ     = -S(end,1,:);          % d(-x)/du  =  - ∂x/∂u
gradJ     = gradJ(:);             % columna
end
%% =====================================================================
%               RESTRICCIONES + GRADIENTES ANALÍTICOS
% ======================================================================
function [c,ceq,gradc,gradceq] = con_with_grad(u,y0,pars)
[Y,S]  = simulate_state_sens(u,y0,pars);

% --- Desigualdades:  y(t) >= 0  →  -y <= 0
c      = -Y(:,2);                       % tamaño N×1
St     = squeeze(S(:,2,:));             % (N × nu)
gradc  = -St.';                         % (nu × N)

% --- Igualdad final: y(tf) = 0
ceq        = Y(end,2);                  % escalar
gradceq    = squeeze(S(end,2,:));       % (nu × 1)
end
%% =====================================================================
%        SIMULACIÓN DEL SISTEMA + SENSIBILIDADES   (con cache)
% ======================================================================
function [Y,S] = simulate_state_sens(u,y0,pars)
%   Devuelve:
%       Y  : (N × 4)  estados en los nodos
%       S  : (N × 4 × nu)  sensibilidades
persistent u_last  Y_last  S_last
if ~isempty(u_last) && isequal(u,u_last)
    Y = Y_last;   S = S_last;   return           % reutiliza integración
end

tf      = u(1);
alpha_n = u(2:end);                 % N valores
nu      = numel(u);

z0          = zeros(4 + 4*nu,1);    % estado + sens.
z0(1:4)     = y0;

t_span  = linspace(0,tf,pars.N);    % nodos reales (N tiempos)

odeopts = odeset('RelTol',1e-6,'AbsTol',1e-8);

odefun  = @(t,z) aug_dynamics(t,z,tf,alpha_n,pars);

[~,Z]   = ode45(odefun,t_span,z0,odeopts);

Y           = Z(:,1:4);                           %   N × 4
Sens        = Z(:,5:end);                         %   N × 4·nu
S           = reshape(Sens,[],4,nu);              %   N × 4 × nu

%  Guarda en cache
u_last = u;   Y_last = Y;   S_last = S;
end
%% =====================================================================
%                 DINÁMICA AUMENTADA  (estado+sens.)
% ======================================================================
function dzdt = aug_dynamics(t,z,tf,alpha_n,pars)
% --- desempaqueta --------------------------------------------
N       = pars.N;
nu      = N + 1;
y       = z(1:4);                      % estado
S       = reshape(z(5:end),4,nu);      % 4×nu sensibilidades

x      = y(1);    y_pos = y(2);        %#ok<NASGU>
v      = y(3);    gamma  = y(4);

% --- α(t), dα/dtf y pesos dα/dα_k -----------------------------
[alpha, w, dalpha_dtf] = alpha_info(t,tf,alpha_n,pars);

% ---- COEFICIENTES AERODINÁMICOS ------------------------------
g   = pars.g;       rho = pars.rho;   Sw = pars.Sw;
m   = pars.m;       Cd0 = pars.Cd0;   k  = pars.k;
Cl0 = pars.Cl0;     Clalpha = pars.Clalpha;

Cl      = Cl0 + Clalpha*alpha;
Cd      = Cd0 + k*Cl.^2;

D       = 0.5*rho*Sw*v^2*Cd;
L       = 0.5*rho*Sw*v^2*Cl;

% --- DINÁMICA f(y,α) ------------------------------------------
f = [ ...
      v*cos(gamma); ...
      v*sin(gamma); ...
     -g*sin(gamma)  - D/m; ...
     -(g/v)*cos(gamma) + L/(m*v) ];

%% ---------------   MATRIZ A = ∂f/∂y   -------------------------
A = zeros(4,4);
A(1,3) =  cos(gamma);
A(1,4) = -v*sin(gamma);

A(2,3) =  sin(gamma);
A(2,4) =  v*cos(gamma);

c1      = 0.5*rho*Sw*Cd;
A(3,3)  = -(1/m) * 2*c1*v;
A(3,4)  = -g*cos(gamma);

k2      = (0.5*rho*Sw*Cl)/m;
A(4,3)  =  g*cos(gamma)/v^2  + k2;
A(4,4)  =  (g/v)*sin(gamma);

%% -------- B = ∂f/∂u_j   (α_k  y  t_f) -------------------------
%  d f / d α
dCd_dalp  = 2*k*Cl*Clalpha;
dD_dalp   = 0.5*rho*Sw*v^2*dCd_dalp;
dL_dalp   = 0.5*rho*Sw*v^2*Clalpha;

df_dalp = [ 0; 0;
           -(1/m)*dD_dalp;
            (1/(m*v))*dL_dalp ];

B = zeros(4,nu);
%  primera columna → t_f
B(:,1)   = df_dalp * dalpha_dtf;     % 4×1
%  columnas 2..N+1 → α_k
B(:,2:end) = df_dalp .* w;           % 4×N   (broadcast)

%% ------------  Ecuación de sensibilidades --------------------
dSdt = A*S + B;                      % 4×nu

%% ------------  Empaqueta derivadas ---------------------------
dzdt = [ f ;  dSdt(:) ];
end
%% =====================================================================
%        α(t), pesos dα/dα_k y derivada dα/dt_f  (malla s = t/tf)
% ======================================================================
function [alpha, w, dalpha_dtf] = alpha_info(t,tf,alpha_n,pars)
N   = pars.N;
s   = min(max(t/tf,0),1);            % s ∈ [0,1]
% nodos uniformes en s
i     = floor(s*(N-1))+1;            % índice base (1…N)
if i>=N
    i       = N-1;   lambda = 1;     % extremo derecho
elseif i<=0
    i       = 1;     lambda = 0;     % extremo izquierdo
else
    s_i     = (i-1)/(N-1);
    lambda  = (s - s_i) * (N-1);     % ∈ [0,1]
end
% pesos lineales
w           = zeros(1,N);
w(i)        = 1 - lambda;
w(i+1)      = lambda;

alpha       = w*alpha_n.';           % α(t)

% ←—— derivada respecto tf ———————————————
dalpha_ds   = (N-1)*(alpha_n(i+1) - alpha_n(i));   % pendiente local
dalpha_dtf  = - (s/tf) * dalpha_ds;                %  dα/dtf
end
