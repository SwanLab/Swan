function [f,gradf] = objfungrad(x,num_cons,phisical_cons)
h = num_cons;
f = -x(4.*(h-1)+1);
% Gradient of the objective function:
gradf = zeros(4.*h+2,1);
gradf(4.*(h-1)+1,1) = -1;

end
