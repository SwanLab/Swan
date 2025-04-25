function main_constraint_nodrag

persistent u_history J_history;
clc;  close all

% Physical constants

g = 9.81; % m/s^2
lb = [0 0.01]; % [m/s rad]
ub = [10 pi/2]; % [m/s rad]

% Initial conditions

x1_0 = 0; % m
x2_0 = 0; % m
v0 = 15; % m/s
gamma0 = deg2rad(40); % degrees
t0 = 0; % s
tf = 2; % s
u0 = [tf;gamma0];


cost = @(u) f_cost(u,g,t0,v0,x1_0,x2_0);
constraint = @(u) groundConstraint(u,g,t0,v0,x1_0,x2_0);

%Check gradients
valid = checkGradients(cost,u0)
valid = checkGradients(constraint,u0)

% Fmincon
options = optimoptions("fmincon","OutputFcn",@store_fmincon,"Algorithm","sqp",'DerivativeCheck','on');
[u_opt, J_opt, exitflag] = fmincon(cost,u0, ...
    [],[],[],[],lb,ub,@(u) nonlcon(u,g,t0,v0,x1_0,x2_0),options);


% Optimal solution

y0 = [x1_0 x2_0 v0 u_opt(2)];
t_span = linspace(t0,u_opt(1),100);

[~,y] = ode45(@(t,y) dynamics(y,g), t_span, y0);

x1_opt = y(:,1);
x2_opt = y(:,2);
tf_opt = u_opt(1);
gamma0_opt = u_opt(2);

% Postprocess

disp(["Maximum distance [m] = ",num2str(x1_opt(end))])
disp(["Initial angle [°]: ",num2str(rad2deg(gamma0_opt))])
disp(["Final Time [s]: ",num2str(tf_opt)])

% Trajectory

figure
plot(x1_opt,x2_opt,'b-','LineWidth',2)
xlabel("Horizontal distance [m]")
ylabel("Vertical distance [m]")
grid on

% Evolution of u

figure;
plot(u_history(:,1), 'b-o', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Initial velocity u [m/s]');
grid on

figure
plot(u_history(:,2)*180/pi, 'r-o', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Initial angle \gamma ');
grid on;

% Evolution of cost function

figure;
plot(J_history,'o-','LineWidth',2)
xlabel('Iteration')
ylabel('Cost function J')
grid on

% Save results

vars = whos; 

nodrag_results = struct();
for i = 1:length(vars)
    var_name = vars(i).name;  
    nodrag_results.(var_name) = eval(var_name); 
end


function stop = store_fmincon(u,optimValues,state)



if isempty(u_history)
    u_history = [];
    J_history = [];
end

if optimValues.iteration == 0
    u_history = [];
    J_history = [];
end


switch state
    case 'iter'
        % Save value of u for every iteration
        u_history = [u_history; u(:)'];
        J_history = [J_history; optimValues.fval];
end

% Saving in a global variable
assignin('base','u_iterations',u_history);
assignin('base','J_iterations',J_history);

stop = false; % Don't stop the optimization

end


save('nodrag_data.mat', 'nodrag_results');

end


function [c,ceq,Dc,Dceq] = nonlcon(u,g,t0,v0,x1_0,x2_0)
c = [];
Dc = [];
[ceq,Dceq] = groundConstraint(u,g,t0,v0,x1_0,x2_0);
end


function [J,gradJ] = f_cost(u,g,t0,v0,x1_0,x2_0)

y0 = [x1_0 x2_0 v0 u(2)];
t_span = linspace(t0,u(1),1000);

[~,y] = ode45(@(t,y) dynamics(y,g), t_span, y0);

J = -y(end,1);

% compute q

pT = [1 0 0 0];

[~,p] = ode45(@(t,p) p_ode(t,p,y,g,t_span),flip(t_span),pT);

q = p(end,:);

% gradient of cost function with respect to u

dydt_final = dynamics(y(end,:)',g);

gradJ(1,1) = -dydt_final(1);
gradJ(2,1) = -q(1,4);

end


function [ceq,Dceq] = groundConstraint(u,g,t0,v0,x1_0,x2_0)

y0 = [x1_0 x2_0 v0 u(2)]';
t_span = linspace(t0,u(1),1000);
[~,y] = ode45(@(t,y) dynamics(y,g), t_span, y0);

% Constraint y2(T) = 0

ceq = y(end,2);


% Derivatives

pT = [0 -1 0 0];

[~,p] = ode45(@(t,p) p_ode(t,p,y,g,t_span),flip(t_span),pT);

q = p(end,:);

dydt_final = dynamics(y(end,:)',g);

Dceq(1,1) = dydt_final(2);
Dceq(2,1) = -q(4);



end


function dpdt = p_ode(t,p,y,g,t_span)

v = interp1(t_span, y(:,3), t);  % velocity at time t
gamma = interp1(t_span, y(:,4), t);

J = [0 0 0 0;
     0 0 0 0;
     cos(gamma) sin(gamma) 0 (g/v^2)*cos(gamma);
     -v*sin(gamma) v*cos(gamma) -g*cos(gamma) (g/v)*sin(gamma)];

dpdt = -J*p;


end


function dydt = dynamics(y,g)

% State variables
v = y(3);
gamma = y(4);

% Dynamics equation

dydt = zeros(4,1);
dydt(1) = v*cos(gamma);
dydt(2) = v*sin(gamma);
dydt(3) = -g*sin(gamma);
dydt(4) = -(g/v)*cos(gamma);

end