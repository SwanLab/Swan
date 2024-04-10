function [f,gradf] = objfungrad(x,h)
h = 300;
f = -x(4.*h+1);
% Gradient of the objective function:
gradf = zeros(4.*h+2,1);
gradf(4.*h+1,1) = -1;

end
