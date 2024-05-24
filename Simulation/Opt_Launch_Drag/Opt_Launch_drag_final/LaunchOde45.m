function [t,y] = LaunchOde45(gamma0,tf,num_cons,phisical_cons)
g = phisical_cons(1,1);
v_0 =phisical_cons(1,2);
x_1_0 =phisical_cons(1,3);
x_2_0 =phisical_cons(1,4);
t_0 =phisical_cons(1,5);

tspan = [t_0 tf];
y0 = [x_1_0 x_2_0 v_0 gamma0];

%Ode calculation:
[t,y] = ode45(@(t,y) odefcn(t,y,num_cons,phisical_cons),tspan,y0);




end

