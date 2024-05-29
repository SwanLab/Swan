function D_uL = Launch_DuF_calc(T_f,p,tp,y,ty,PI,ALPHA,num_cons,Initial_cons,Engine_cons,Aircraft_cons,Aero_cons,Earth_cons)

Solution_size = num_cons(1,3);
t_0 = Initial_cons(1,6);
t = linspace(t_0,T_f,Solution_size);
y_new = zeros(Solution_size,5);
p_new = zeros(Solution_size,5);

for i = 1:Solution_size
[y_new(i,1),y_new(i,2),y_new(i,3),y_new(i,4),y_new(i,5)] = Interpolation_fun_p(y,ty,t(i));

[p_new(i,1),p_new(i,2),p_new(i,3),p_new(i,4),p_new(i,5)] = Interpolation_fun_p(p,tp,t(i));
end

D_uL = zeros(Solution_size,2);

for i = 1:Solution_size
x1 = y_new(i,1);
x2 = y_new(i,2);
v = y_new(i,3);
gamma = y_new(i,4);
weight = y_new(i,5);

PI_Ajusted = PI(1,i);
ALPHA_Ajusted = ALPHA(1,i);

p_vector = p_new(i,:);

%% Previous calculations:
%Temperature:
Temp = Temperature_fun(x2,Earth_cons);

%AirDensity:
rho = AirDensity_fun(x2,Earth_cons);

%Mach number:
M = Mach_fun(Temp,v,Earth_cons);

%Cl_alpha coefficent:
CL_alpha = CL_alpha_fun(M,Aircraft_cons,Aero_cons);

%CD_0 drag coefficent:
CD_0 = CD_0_fun(M,Aero_cons);

%Drag coefficent:
C_D = C_D_fun(CL_alpha,CD_0,ALPHA_Ajusted,Aero_cons);

%Lift coefficent:
C_L = CL_alpha.*ALPHA_Ajusted;

%Air density before entering the engine:
rho_eng = rho_eng_fun(M,rho,Earth_cons);

%Thrust calculation:
Thrust = Thrust_fun(rho_eng,PI_Ajusted,Engine_cons,Earth_cons);

%Lift force:
Lift = Lift_fun(C_L,rho,v,Aircraft_cons);

%Drag force:
Drag = Lift_fun(C_D,rho,v,Aircraft_cons);

%% Calculation previous derivatives:
Thrust_D_PI = Thrust_fun_D_PI(rho_eng,PI_Ajusted,Engine_cons,Earth_cons);

Lift_D_ALPHA = Lift_fun_D_ALPHA(v,rho,ALPHA_Ajusted,CL_alpha,Aircraft_cons);

Drag_D_ALPHA = Drag_fun_D_ALPHA(v,rho,ALPHA_Ajusted,CL_alpha,Aircraft_cons,Aero_cons);

%% Matrix Fu calculation:
% Constants:
g = Earth_cons(1,1);
C_T = Engine_cons(1,3);
Du_F = zeros(5,2);

% First and second term = 0

% Third term (v):
Du_F(3,1) = (g./weight).*Thrust_D_PI;

Du_F(3,2) = -(g./weight).*Drag_D_ALPHA;

% Fourth term (gamma):
Du_F(4,2) = (g./(v.*weight)).*Lift_D_ALPHA;

% Fifth term (weight):
Du_F(5,1) = - C_T.*Thrust_D_PI;


%% Matrix Du_L calculation:
D_uL(i,:) = p_vector*Du_F;
end


%% Interpolation ODE
%Correction for interpolation of values:




end

