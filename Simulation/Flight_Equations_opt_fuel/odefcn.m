function dydt = odefcn(t,y,PI,ALPHA,T_f,num_cons,Initial_cons,Engine_cons,Aircraft_cons,Aero_cons,Earth_cons)
%% Variables asignation:
x1 = y(1);
x2 = y(2);
v = y(3);
gamma = y(4);
weight = y(5);
%Solution ajustment:
[ALPHA_Ajusted, PI_Ajusted] = Interpolation_fun(ALPHA,PI,T_f,t,num_cons);

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

%% Constants:
g = Earth_cons(1,1);
C_T = Engine_cons(1,3);


dydt = zeros(5,1);

dydt(1) = v.*cos(gamma);
dydt(2) = v.*sin(gamma);
dydt(3) = g.*((Thrust-Drag)./weight-sin(gamma));
dydt(4) = (g./v).*((Lift./weight)-cos(gamma));
dydt(5) = -C_T.*Thrust;

end


