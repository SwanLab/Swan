clc
clear all

% Define the objective function


%Phisical Constants:
g = 9.81;
x_1_0 = 0;
x_2_0 = 0;
t_0 = 0;
t_f = 10;
%Numerical constant:
h = 10000;
weight = 3;


%Phisical data of projectail:
Air_density = 1.225; %kg/m^3;
Surface = 0.01889; %m^2;
Drag_coeffiecent = 0.35; %no units.
Mass = 45; %kg

K_projectile = 0.5.*Air_density.*Surface.*Drag_coeffiecent;

%Numerica assignation:
num_cons = zeros(1,3);
num_cons(1,1) = h;

num_cons(1,3) = weight;

%Phisical asignation:
phisical_cons = zeros(1,9);
phisical_cons(1,1) = g;
phisical_cons(1,2) = 0;
phisical_cons(1,3) = x_1_0;
phisical_cons(1,4) = x_2_0;
phisical_cons(1,5) = t_0;
phisical_cons(1,6) = K_projectile;
phisical_cons(1,7) = Mass;
phisical_cons(1,8) = weight;
phisical_cons(1,9) = t_f;


tolerance = 0.5;
x = zeros(2,1);
x(1,1) = 0.8;
x(2,1) = 100;
previous_fun = 0;

% Initial step sizes
alpha_gamma0 = 0.01;  % Step size for gamma0
alpha_v0 = 0.05;      % Step size for v0

maxIter = 3000;
% Assume 'grad_gamma0' and 'grad_v0' are your calculated gradients for gamma0 and v0
% and 'gamma0' and 'v0' are your current values of these parameters.
Convergence_rate = zeros(maxIter);
for iter= 1:maxIter


    %bounds
    lb = [0 0.01];
    
    ub = [0.5.*pi;inf];

    %Function:
    [fun, grad] = objfungrad(x, num_cons, phisical_cons);  % Compute objective and gradient
    fprintf('Current objective value: %f\n', fun);
    fprintf('Gradient: %f %f\n', grad(1), grad(2));
    gamma0= x(1);
    v0 = x(2);
    % Update parameters
    gamma0 = gamma0 - alpha_gamma0 * grad(1) / (1 + abs(grad(1)));
    v0 = v0 - alpha_v0 * grad(2)/(1 + abs(grad(2)));

    % Optional: Adjust step size dynamically based on some criteria
    if mod(iter, 10) == 0  % Adjust step size every 10 iterations
        alpha_gamma0 = alpha_gamma0 * 0.9;  % Reduce step size by 10%
        alpha_v0 = alpha_v0 * 0.9;
    end

    % Check for convergence (this is a simple placeholder condition)
    if abs(fun - previous_fun) < tolerance
        last_iter = iter;
        break;
    end
    previous_fun = fun;
    x(1) = gamma0;
    x(2) = v0;
    for i = 1:2
        if x(i)>ub(i)
            x(i) = ub(i);
        elseif x(i)<lb(i)
            x(i) = lb(i);
        end
    end
    fun_ant = fun;
    Convergence_rate(iter) = fun;
end



%options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
% Solve the optimization problem
%[x, fval, exitflag, output]= fmincon(fun,x0,[],[],[],[],lb,ub,[],options);

%Representation:
gamma0 = x(1);
tf = x(2);


[t,y]=LaunchOde45(gamma0,tf,num_cons,phisical_cons);


Data_representation(y,t,phisical_cons,Convergence_rate,last_iter)
disp(['Number of iterations: ', num2str(iter)]);
