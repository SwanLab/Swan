%Constants:
g = 9.81;
v_0 = 10;
x_1_0 = 0;
x_2_0 = 0;
t_0 = 0;
h = 300;
% Define the objective function
fun = @(x) objfungrad(x,h);
nonlincon =@nlcon(x,h);




num_cons = zeros(1,1);
num_cons(1,1) = h;
phisical_cons = zeros(1,5);
phisical_cons(1,1) = g;
phisical_cons(1,2) = v_0;
phisical_cons(1,3) = x_1_0;
phisical_cons(1,4) = x_2_0;
phisical_cons(1,5) = t_0;


% Initial guess
x0 = zeros(1, 4.*h+2);
x0(1,4.*h+2) = 1.2;
options = optimoptions('fmincon','GradObj','on','GradConstr','on');




%bounds
lb = zeros(4.*h+2,1);

for i=1:h
    lb(4.*i,1) = -inf;
end


% Solve the optimization problem

x = fmincon(fun,x0,[],[],[],[],lb,[],nonlincon,options)





