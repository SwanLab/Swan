function [f,gradf] = objfungrad(x,num_cons,phisical_cons)

%Numerical values:
alpha = num_cons(1,2);
Weight = num_cons(1,3);

% Variables:
gamma0 = x(1);
v_0 = x(2);

%% Max distance:
[t,y]=LaunchOde45(gamma0,v_0,num_cons,phisical_cons);

f = -y(end,1)+0.5.*Weight.*(y(end,2)).^2;
%% Calculation of p value 
ty = t;
[t,p]=LaunchOde45_P_value(gamma0,v_0,num_cons,phisical_cons,y,ty);


q = p(end,:);
D_uL = -q(1,4);
D_uL_2 = -q(1,3);

%% Gradient of the objective function:
gradf = [D_uL,D_uL_2];
% For debugging:
disp(['Current state: ', num2str(y(end,:))]);
disp(['Current objective value: ', num2str(f)]);
disp(['Gradient: ', num2str(gradf)]);

end

