function [f,gradf] = objfungrad(x,num_cons,phisical_cons)

gamma0 = x(1);
tf = x(2);


%Numerical constants:
Weight = num_cons(1,2);
alpha = num_cons(1,3);
%% Max time:


%% Max distance:
[t,y]=LaunchOde45(gamma0,tf,num_cons,phisical_cons);

%Optimization function:
f = -y(end,1)+0.5.*Weight.*(y(end,2)).^2;
%% Calculation of p value 
%Transformation of previous time
ty = t;
[t,p]=LaunchOde45_P_value(gamma0,tf,num_cons,phisical_cons,y,ty);

% q calculation:
q = p(1,:);
D_uL = -q(1,4).*alpha;

%% Gradient of the objective function:
gradf = [D_uL;0];



end

