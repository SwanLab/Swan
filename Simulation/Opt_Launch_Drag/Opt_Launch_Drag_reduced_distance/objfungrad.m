function f = objfungrad(x,num_cons,phisical_cons)
h = num_cons;

gamma0 = x(1);

tf = x(2);
y = launch(gamma0,tf,num_cons,phisical_cons);

f = -y(1,end);
%f = -tf;

% Gradient of the objective function:
%gradf = [0;0];
end

