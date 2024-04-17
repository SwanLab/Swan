function [f,gradf] = objfungrad(x,num_cons,phisical_cons)
h = num_cons;
f = -x(2);
% Gradient of the objective function:
gradf = [0;-1];


end

