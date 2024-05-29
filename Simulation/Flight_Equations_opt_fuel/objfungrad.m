function [f,gradf] = objfungrad(x,num_cons,Initial_cons,Engine_cons,Aircraft_cons,Aero_cons,Earth_cons)

PI = x(1,:);
ALPHA = x(2,:);
T_f = x(3,1);

%Numerical constants:
weight_numerical = num_cons(1,1);
alpha_numerical = num_cons(1,2);
h_objective = num_cons(1,4);
Solution_size = num_cons(1,3);

%% Max time:


%% Max distance:
[t,y]=LaunchOde45(PI,ALPHA,T_f,num_cons,Initial_cons,Engine_cons,Aircraft_cons,Aero_cons,Earth_cons);

%Optimization function:
f = -y(end,5)+0.5.*weight_numerical.*(y(end,2)-h_objective).^2;


%% Calculation of p value 
%Transformation of previous time
ty = t;
[t,p]=LaunchOde45_P_value(PI,ALPHA,T_f,num_cons,Initial_cons,Engine_cons,Aircraft_cons,Aero_cons,Earth_cons,y,ty);

tp = t;

D_uL = Launch_DuF_calc(T_f,p,tp,y,ty,PI,ALPHA,num_cons,Initial_cons,Engine_cons,Aircraft_cons,Aero_cons,Earth_cons);


%% Gradient of the objective function:
gradf = zeros(Solution_size,3);
gradf(:,1:2) = D_uL;

gradf = alpha_numerical.*gradf;

gradf = transpose(gradf);
 



end

