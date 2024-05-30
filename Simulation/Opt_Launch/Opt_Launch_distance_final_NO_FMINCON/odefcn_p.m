function dpdt = odefcn_p(t,p,num_cons,phisical_cons,y,ty)
g = phisical_cons(1,1);


x_1 = interp1(ty,y(:,1),t);
x_2 = interp1(ty,y(:,2),t);
v_ant = interp1(ty,y(:,3),t);
gamma_ant = interp1(ty,y(:,4),t);

J = [0 0 0 0;
     0 0 0 0;
     cos(gamma_ant) sin(gamma_ant) 0 (g./(v_ant.^2)).*cos(gamma_ant);
     -v_ant.*sin(gamma_ant) v_ant.*cos(gamma_ant) -g.*cos(gamma_ant) (g./v_ant).*sin(gamma_ant)];


dpdt = -J*p;

end


