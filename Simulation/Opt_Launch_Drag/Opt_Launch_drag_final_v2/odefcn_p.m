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

%for i = 1:st
%x_1 = y(i,1);
%x_2 = y(i,2);
%v_ant = y(i,3);
%gamma_ant = y(i,4);


%F(:,:,i) = [0 0 0 0;
 %    0 0 0 0;
 %    cos(gamma_ant) sin(gamma_ant) -2.*(K./m).*v_ant (g./(v_ant.^2)).*cos(gamma_ant);
 %    -v_ant.*sin(gamma_ant) v_ant.*cos(gamma_ant) -g.*cos(gamma_ant) (g./gamma_ant).*sin(gamma_ant)];
%F(:,:,i) = transpose(F(:,:,i));
%end
%J = zeros(4,4);
%X_search = zeros(st,1);
%for i = 1:4
%    for j = 1:4
%        X_search =  F(i,j,:);
%        X_search = squeeze(X_search);
        %X_search = transpose(X_search);
        %ty = transpose(ty);
%        J(i,j) = interp1(ty,X_search,t);
 %   end
%end


%F = interp1(F,ty,t);

dpdt =-J*p;

end


