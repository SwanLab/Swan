clc

% Define the objective function


%Constants:
g = 9.81;
v_0 = 10;
x_1_0 = 0;
x_2_0 = 0;
t_0 = 0;
h = 500;

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
%x0(1,4.*h+2) = 0.5;
%options = optimoptions('fmincon','GradObj');




%bounds
lb = zeros(2,1);
ub = [inf;0.5.*pi];

% Solve the optimization problem

x = fmincon(fun,x0,[],[],[],[],lb,ub,nonlincon);

y = launch(x(2),x(1),num_cons,phisical_cons);


plot(y(1,:),y(2,:))
disp('time')
x(1)
disp('angle')
x(2)