clc
clear all


%% Constants:
g = 9.81;
v_0 = 10;
x_1_0 = 0;
x_2_0 = 0;
t_0 = 0;

%Numerical constants:
h = 10000;
weight = 100;
alpha = 0.5;

%Phisical data of projectail:
Air_density = 1.225; %kg/m^3;
Surface = 2; %m^2;
Drag_coeffiecent = 0.4; %no units.
Mass = 100; %kg

K_projectile = 0.5.*Air_density.*Surface.*Drag_coeffiecent;

%% Value assignation:

num_cons = zeros(1,3);
num_cons(1,1) = h;
num_cons(1,2) = weight;
num_cons(1,3) = alpha;
phisical_cons = zeros(1,7);
phisical_cons(1,1) = g;
phisical_cons(1,2) = v_0;
phisical_cons(1,3) = x_1_0;
phisical_cons(1,4) = x_2_0;
phisical_cons(1,5) = t_0;
phisical_cons(1,6) = K_projectile;
phisical_cons(1,7) = Mass;


fun = @(x) objfungrad(x,num_cons,phisical_cons);


% Initial guess
x0 = zeros(2,1);
x0(1,1) = 1;


%bounds
lb = zeros(2,1);
ub = [0.5.*pi;inf];

% Solve the optimization problem
[x, fval, exitflag, output]= fmincon(fun,x0,[],[],[],[],lb,ub);

%Representation:
gamma0 = x(1);
tf = x(2);


[t,y]=LaunchOde45(gamma0,tf,num_cons,phisical_cons);

Data_representation(y,t,phisical_cons,output)
