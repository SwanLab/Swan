clc

% Define the objective function


%Constants:
g = 9.81;
v_0 = 10;
x_1_0 = 0;
x_2_0 = 0;
t_0 = 0;

%Phisical data of projectail:
Air_density = 1.225; %kg/m^3;
Surface = 2; %m^2;
Drag_coeffiecent = 0.4; %no units.
Mass = 5; %kg

K_projectile = 0.5.*Air_density.*Surface.*Drag_coeffiecent;


h = 500;

%num_cons = zeros(1,1);
num_cons = h;
phisical_cons = zeros(1,5);
phisical_cons(1,1) = g;
phisical_cons(1,2) = v_0;
phisical_cons(1,3) = x_1_0;
phisical_cons(1,4) = x_2_0;
phisical_cons(1,5) = t_0;
phisical_cons(1,6) = K_projectile;
phisical_cons(1,7) = Mass;

fun = @(x) objfungrad(x,num_cons,phisical_cons);
nonlincon = @(x) nlcon(x,num_cons,phisical_cons);

% Initial guess
x0 = zeros(2,1);
x0(1,1) = 1;
%x0(1,4.*h+2) = 0.5;
%options = optimoptions('fmincon','GradObj');




%bounds
lb = zeros(2,1);
ub = [0.5.*pi;inf];

% Solve the optimization problem

x = fmincon(fun,x0,[],[],[],[],lb,ub,nonlincon);

y = launch(x(1),x(2),num_cons,phisical_cons);

time = zeros(1,h);
time(1,1) = t_0;

for i = 2:h
time(1,i) = time(1,i-1) + x(2)./h;

end

figure(1)
plot(y(1,:),y(2,:))
title Trajectory
xlabel [m]
ylabel [m]
axis equal
xlim([0 inf])
ylim([0 3])
disp('time')
x(2)
disp('angle')
x(1)


[c,ceq] = nlcon(x,num_cons,phisical_cons)

figure(2)
plot(time(1,:),y(3,:))
title Speed
xlabel [s]
ylabel [m/s]

figure(3)
plot(time(1,:),y(4,:))
title Gamma
xlabel [s]
ylabel [rad]

