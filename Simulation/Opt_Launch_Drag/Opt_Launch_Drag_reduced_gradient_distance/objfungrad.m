function [f,gradf]= objfungrad(x,num_cons,phisical_cons)
h = num_cons;
alpha = 1;
gamma0 = x(1);
tf = x(2);
g = phisical_cons(1,1);
t_0=phisical_cons(1,5);
K_projectile = phisical_cons(1,6);
Mass = phisical_cons(1,7);



y = launch(gamma0,tf,num_cons,phisical_cons);

f = -y(1,end);
%f = -tf;
%% Calculation of p value 
[sy, st] = size(y);
p = zeros(sy,st);
p(1,end) = -1;
p(2,end) = -y(2,end);
Inc_time = ((tf-t_0)./h);
for i = st:-1:2 
    x_1 = y(1,i);
    x_2 = y(2,i);
    v_ant = y(3,i);
    gamma_ant = y(4,i);

    F_matrix = [0 0 0 0;
                0 0 0 0;
                cos(gamma_ant) sin(gamma_ant) -2.*(K_projectile./Mass).*v_ant (g./(v_ant.^2)).*cos(gamma_ant);
                -v_ant.*sin(gamma_ant) v_ant.*cos(gamma_ant) -g.*cos(gamma_ant) (g./gamma_ant).*sin(gamma_ant)];
    p_dot = -F_matrix*p(:,i);
    
    p(:,i-1) = p(:,i)-p_dot.*Inc_time;

end

q = p(:,1);
D_uL = -q(4,1).*alpha;

%% Gradient of the objective function:
gradf = [D_uL;0];


% Gradient of the objective function:
%gradf = [0;0];
end

