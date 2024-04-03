function [c, ceq,gradc,gradceq] = nlincon2(x)
gamma0 = x(1);
tf = x(2);
y = launch(y0,gamma0,tf);
ceq =y(2,end) - 0;
c = [];
gradc = [];
gradceq = [];


end