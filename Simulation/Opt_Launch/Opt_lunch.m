% Define the objective function
fun = @objfungrad;

%Constants:
g = 9.81;
v_0 = 10;
x_1_0 = 0;
x_2_0 = 0;
t_0 = 0;
h = 3;

% Initial guess
x0 = zeros(1, 14); %
options = optimoptions('fmincon','GradObj','on','GradConstr','on');



nonlincon =@nlcon;

%bounds
lb = 0.001 * ones(4);

% Solve the optimization problem

x = fmincon(fun,x0,[],[],[],[],lb,[],nonlincon,options)





