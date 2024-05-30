function dpdt = odefcn_p(t,p,num_cons,phisical_cons,y,ty)
g = phisical_cons(1,1);
K = phisical_cons(1,6);
m = phisical_cons(1,7);


x_1 = interp1(ty,y(:,1),t);
x_2 = interp1(ty,y(:,2),t);
v_ant = interp1(ty,y(:,3),t);
gamma_ant = interp1(ty,y(:,4),t);

J = [0 0 0 0;
     0 0 0 0;
     cos(gamma_ant) sin(gamma_ant) -2.*(K./m).*v_ant (g./(v_ant.^2)).*cos(gamma_ant);
     -v_ant.*sin(gamma_ant) v_ant.*cos(gamma_ant) -g.*cos(gamma_ant) (g./v_ant).*sin(gamma_ant)];


dpdt = -J*p;

end


