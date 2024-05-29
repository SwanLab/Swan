function dpdt = odefcn_p(t,p,PI,ALPHA,T_f,num_cons,Initial_cons,Engine_cons,Aircraft_cons,Aero_cons,Earth_cons,y,ty)
%% Interpolation ODE
%Correction for interpolation of values:
[x1,x2,v,gamma,weight] = Interpolation_fun_p(y,ty,t);

%Correction of variables
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

%% Calculation previous derivatives:
%Derivative air density with height
D_rho_D_h = Derivative_airdensity_d_h(x2,Earth_cons);

%Derivative Temperature with height
D_Temp_D_h = Derivative_Temp_d_h(x2,Earth_cons);

%Derivative of Mach number, height and airspeed.
[D_M_D_h, D_M_D_v] = Derivative_Mach(v,D_Temp_D_h,Temp,Earth_cons);

%Derivative of airdensity before entering the engine, height and airspeed.
[D_rho_eng_D_h, D_rho_eng_D_v] = Derivative_airdensity_eng_d_h(M,rho,Earth_cons,D_M_D_v,D_M_D_h,D_rho_D_h);

%Derivative of CD_0 coefficent, height and airspeed.
[D_C_D_0_D_h, D_C_D_0_D_v] = Derivative_C_D_0_at_mach(M,Aero_cons,D_M_D_v,D_M_D_h);

%Derivative of CL_alpha coefficent, height and airspeed.
[D_CL_alpha_D_h, D_CL_alpha_D_v] = Derivative_CL_alpha_at_mach(M,Aero_cons,Aircraft_cons,D_M_D_v,D_M_D_h);

%Derivative of CL coefficent, height and airspeed.
[D_CL_D_h, D_CL_D_v] = Derivative_CL(D_CL_alpha_D_h,D_CL_alpha_D_v,ALPHA_Ajusted,CL_alpha);

%Derivative of CD coefficent, height and airspeed.
[D_CD_D_h, D_CD_D_v] = Derivative_CD(D_CL_alpha_D_h,D_CL_alpha_D_v,D_C_D_0_D_h,D_C_D_0_D_v,ALPHA_Ajusted,CL_alpha,Aero_cons);

%Derivative of Lift:
[Lift_D_h, Lift_D_v] = Derivative_Lift(v,rho,C_L,Aircraft_cons,D_CL_D_v,D_CL_D_h,D_rho_D_h);

%Derivative of Drag:
[Drag_D_h, Drag_D_v] = Derivative_Lift(v,rho,C_D,Aircraft_cons,D_CD_D_v,D_CD_D_h,D_rho_D_h);

%Derivative of Thrust:
[Thrust_D_h, Thrust_D_v] = Derivative_Thrust(PI_Ajusted,Engine_cons,Earth_cons,D_rho_eng_D_v,D_rho_eng_D_h);

%% Creation of derivation matrix D_Y F(u,y):
%Constants:
g = Earth_cons(1,1);
C_T = Engine_cons(1,3);
Fy = zeros(5,5);
%First term x1:

%Where 1, 2 and 5 = 0
Fy(1,3) = cos(gamma);
Fy(1,4) = - v.*sin(gamma);

%Second term x2 or h:

%Where 1, 2 and 5 = 0
Fy(2,3) = sin(gamma);
Fy(2,4) = v.*cos(gamma);

%Third term v:

%Where 1 = 0
Fy(3,2) = (g./weight).*(Thrust_D_h-Drag_D_h);


Fy(3,3) = (g./weight).*(Thrust_D_v-Drag_D_v);

Fy(3,4) = -g.*cos(gamma);

Fy(3,5) = -(g.*(Thrust-Drag))./(weight.^2);

%Fourth term gamma:

%Where 1 = 0
Fy(4,2) = (g./(v.*weight)).*Lift_D_h;


Fy(4,3) = -(g./(v.^2)).*((Lift./weight)-cos(gamma)) + (g./(v.*weight)).*Lift_D_v;

Fy(4,4) = (g./v).*sin(gamma);

Fy(4,5) = -(g./v).*(Lift./(weight.^2));

%Fifth term weight:

%Where 1, 4 and 5= 0
Fy(5,2) = -C_T.*Thrust_D_h;


Fy(5,3) = -C_T.*Thrust_D_v;




dpdt = -Fy*p;

end


