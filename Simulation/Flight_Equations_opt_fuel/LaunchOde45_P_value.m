function [t,p] = LaunchOde45_P_value(PI,ALPHA,T_f,num_cons,Initial_cons,Engine_cons,Aircraft_cons,Aero_cons,Earth_cons,y,ty)
% Constants:
t_0 = Initial_cons(1,6);
weight_numerical = num_cons(1,1);
h_objective = num_cons(1,4);

tspan = [T_f t_0];
p0 = [0 weight_numerical.*(h_objective- y(2,end)) 0 0 1];

[t,p] = ode45(@(t,p) odefcn_p(t,p,PI,ALPHA,T_f,num_cons,Initial_cons,Engine_cons,Aircraft_cons,Aero_cons,Earth_cons,y,ty),tspan,p0);



end

