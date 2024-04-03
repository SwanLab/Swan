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
options = optimoptions('fmincon','GradObj','on','GradConstr','on','Display','iter');



nonlincon =@nlcon;

%bounds
lb = [0;0;0;-inf;...
    0;0;0;-inf;...
    0;0;0;-inf;...
    0;0];

% Solve the optimization problem

x = fmincon(fun,x0,[],[],[],[],lb,[],nonlincon,options);
%
% 1. %code launch function
% 2. simulate Launch using launch function


% (Full problem)
% 4. nlcon amb loop 

% (Reduced problem)
% 5. Resoldre amb nomes gamma0,tf in Opt_lunch2
% cost [0 -1] i gradient
% nlcon using launch function using y(2,end)
% 6. Write in Matcha Reduced problem