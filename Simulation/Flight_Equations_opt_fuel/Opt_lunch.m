clc
clear all


%% Constants:
%Initial constants
v_0 = 120; % 360 km/h
x_1_0 = 0;
x_2_0 = 0;
gamma_0 = 0;
Takeoff_mass = 78000; %KG % Fuel weight = 15 500 kg  This is the maximum take of weight
t_0 = 0;

Takeoff_weight = Takeoff_mass.*9.81;

% Engine data: IAE V2522 A5
Thrust_max = 111205; % N each engine.
Thrust_idle = Thrust_max.*0.65; % 65% of max thrust
C_T = 9.91720.*10^(-6); %(kg/s)/N
N_engine = 2; %Number of engines

Thrust_max = N_engine.*Thrust_max;
Thrust_idle = N_engine.*Thrust_idle;

%Airbus 320 specifications:
S = 122.6; %m^2 wing surface
Wingspan = 38.8; %m 
AR = ((Wingspan)^2)/S;


%Aerodinamic coefficents:
e =  0.9220;   % Osbald factor
K =  0.0334;   % polar curve factor
C_l_alpha_0 = 2.*pi; % Clalpha of profile
C_D_0.a = 0.023;
C_D_0.b =-0.001533333;
C_D_0.c = 0.020444444;
C_D_0.d =-0.072577778;
C_D_0.e = 0.079733333;

%Earth phisical constants:

g = 9.81; %Gravity m/s2
rho_0 = 1.225; %kg/m3
rho_11 = 0.3636; %kg/m3
L_temp = 0.0065;% K/m
Temperature_0 = 288.15; %K
Temperature_11 = 216.65; %K
r_gas = 287; %J/(kg·K)
gamma_air = 1.4; % diatomic gas



%Numerical constants:

weight_numerical = 10000000;
alpha_numerical = 0.8;
Solution_size = 1000; %A solution every 0.25 seconds with max time
h_objective = 10000; %In ft;
h_objective = 0.3048.*h_objective; % In meters

%% Value assignation:
%Numerical constants
num_cons = zeros(1,4);
num_cons(1,1) = weight_numerical;
num_cons(1,2) = alpha_numerical;
num_cons(1,3) = Solution_size;
num_cons(1,4) = h_objective;

%Initial constants:
Initial_cons = zeros(1,6);
Initial_cons(1,1) = x_1_0;
Initial_cons(1,2) = x_2_0;
Initial_cons(1,3) = v_0;
Initial_cons(1,4) = gamma_0;
Initial_cons(1,5) = Takeoff_weight;
Initial_cons(1,6) = t_0;

%Engine constants:
Engine_cons = zeros(1,3);
Engine_cons(1,1) = Thrust_max;
Engine_cons(1,2) = Thrust_idle;
Engine_cons(1,3) = C_T;

%Aircraft constants:
Aircraft_cons = zeros(1,3);
Aircraft_cons(1,1) = S;
Aircraft_cons(1,2) = Wingspan;
Aircraft_cons(1,3) = AR;

%Aerodynamic coefficents:
Aero_cons = zeros(1,8);
Aero_cons(1,1) = K;
Aero_cons(1,2) = C_l_alpha_0;
Aero_cons(1,3) = e;
Aero_cons(1,4) = C_D_0.a;
Aero_cons(1,5) = C_D_0.b;
Aero_cons(1,6) = C_D_0.c;
Aero_cons(1,7) = C_D_0.d;
Aero_cons(1,8) = C_D_0.e;

%Earth constants:
Earth_cons = zeros(1,8);
Earth_cons(1,1) = g;
Earth_cons(1,2) = rho_0;
Earth_cons(1,3) = rho_11;
Earth_cons(1,4) = Temperature_0;
Earth_cons(1,5) = Temperature_11;
Earth_cons(1,6) = L_temp;
Earth_cons(1,7) = r_gas;
Earth_cons(1,8) = gamma_air;


fun = @(x) objfungrad(x,num_cons,Initial_cons,Engine_cons,Aircraft_cons,Aero_cons,Earth_cons);


% Initial guess

x0 = zeros(3,Solution_size);
x0(1,:) = 1;
x0(2,:) = 0.05.*pi; % 9 degreas of alpha
x0(3,:) = 1000;
%% bounds
%Lower
lb = zeros(3,Solution_size);

lb(2,:)=-(0.5).*pi; 

%Upper
ub = zeros(3,Solution_size);
ub(1,:) = 1;
ub(2,:) = 0.1.*pi;
ub(3,1) = 1800;

options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
% Solve the optimization problem
[x, fval, exitflag, output]= fmincon(fun,x0,[],[],[],[],lb,ub,[],options);

%Representation:
gamma0 = x(1);
tf = x(2);


[t,y]=LaunchOde45(gamma0,tf,num_cons,phisical_cons);

Data_representation(y,t,phisical_cons,output)
