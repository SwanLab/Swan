function [t,y] = LaunchOde45(PI,ALPHA,T_f,num_cons,Initial_cons,Engine_cons,Aircraft_cons,Aero_cons,Earth_cons);
%% Constants:
%Initial_cons
x_1_0 = Initial_cons(1,1);
x_2_0 = Initial_cons(1,2);
v_0 = Initial_cons(1,3);
gamma_0 = Initial_cons(1,4);
W_0 = Initial_cons(1,5);
t_0 = Initial_cons(1,6);

%% Initial values
tspan = [t_0 T_f];
y0 = [x_1_0 x_2_0 v_0 gamma_0 W_0];

%% Ode calculation:
[t,y] = ode45(@(t,y) odefcn(t,y,PI,ALPHA,T_f,num_cons,Initial_cons,Engine_cons,Aircraft_cons,Aero_cons,Earth_cons),tspan,y0);




end

