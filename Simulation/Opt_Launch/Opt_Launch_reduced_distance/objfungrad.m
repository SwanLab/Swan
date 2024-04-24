function [f,gradf] = objfungrad(x,num_cons,phisical_cons)
h = num_cons;
gamma0 = x(1);
tf = x(2);
%% Max time:
%f = -tf;
%Gradient
%gradf = [0;-1];

%% Max distance:
y = launch(gamma0,tf,num_cons,phisical_cons);
f = -y(1,end);
%% Gradient of the objective function:
gradf = [0;0];



end

