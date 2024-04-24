clc
clear all

% Define the objective function


%Constants:
g = 9.81;
v_0 = 10;
x_1_0 = 0;
x_2_0 = 0;
t_0 = 0;
h = 10000;

%num_cons = zeros(1,1);
num_cons = h;
phisical_cons = zeros(1,5);
phisical_cons(1,1) = g;
phisical_cons(1,2) = v_0;
phisical_cons(1,3) = x_1_0;
phisical_cons(1,4) = x_2_0;
phisical_cons(1,5) = t_0;

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
gamma0 = x(1);
tf = x(2);
y = launch(gamma0,tf,num_cons,phisical_cons);

time = zeros(1,h);
time(1,1) = t_0;

for i = 2:h
time(1,i) = time(1,i-1) + tf./h;

end

figure(1)
plot(y(1,:),y(2,:))
title Trajectory
xlabel [m]
ylabel [m]
axis equal
disp('time')
x(2)
disp('angle')
ylim([0 4])
xlim([0 inf])
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

