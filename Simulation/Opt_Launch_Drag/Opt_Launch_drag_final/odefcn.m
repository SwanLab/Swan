function dydt = odefcn(t,y,num_cons,phisical_cons)
g = phisical_cons(1,1);
K = phisical_cons(1,6);
m = phisical_cons(1,7);

x1 = y(1);
x2 = y(2);
v = y(3);
gamma = y(4);
dydt = zeros(4,1);

dydt(1) = v.*cos(gamma);
dydt(2) = v.*sin(gamma);
dydt(3) = -g.*sin(gamma) -(K./m).*v.^2;
dydt(4) = -(g./v).*cos(gamma);

end


