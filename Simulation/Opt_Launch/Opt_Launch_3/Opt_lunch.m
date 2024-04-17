clc

% Define the objective function


%Constants:
g = 9.81;
v_0 = 10;
x_1_0 = 0;
x_2_0 = 0;
t_0 = 0;
h = 100;

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
x0 = zeros(1, 4.*h+2);
x0(1,4.*h+2) = 0.5;
options = optimoptions('fmincon','GradObj','on','GradConstr','on');




%bounds
lb = zeros(4.*h+2,1);

for i=1:(h-1)
    lb(4.*(i-1)+4,1) = -inf;
end


% Solve the optimization problem

x = fmincon(fun,x0,[],[],[],[],lb,[],nonlincon,options);

x1 = zeros(h);
x2 = zeros(h);
v = zeros(h);
gamma = zeros(h);
for i=0:(h-1)
x1(i+1) = x(i+1);
x2(i+1) = x(i+2);
v(i+1) = x(i+3);
gamma(i+1) = x(i+4);
end

plot(x1,x2)




