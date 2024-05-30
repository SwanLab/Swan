function [t,p] = LaunchOde45_P_value(gamma0,v_0,num_cons,phisical_cons,y,ty)
g = phisical_cons(1,1);
v_0 =phisical_cons(1,2);
x_1_0 =phisical_cons(1,3);
x_2_0 =phisical_cons(1,4);
t_0 =phisical_cons(1,5);
tf = phisical_cons(1,9);
weight = phisical_cons(1,8);

tspan = [tf t_0];
p0 = [1 -weight.*y(end,2) 0 0];


options = odeset('RelTol',1e-8,'AbsTol',1e-10);  % More stringent tolerances
%launch_function = odefcn_p(t,y,num_cons,phisical_cons,y,ty);

[t,p] = ode45(@(t,p) odefcn_p(t,p,num_cons,phisical_cons,y,ty),tspan,p0,options);


end

