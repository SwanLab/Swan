function [x1,x2,v,gamma,weight] = Interpolation_fun_p(y,ty,t)
%Calculations:
x1 = interp1(ty,y(:,1),t);
x2 = interp1(ty,y(:,2),t);
v = interp1(ty,y(:,3),t);
gamma = interp1(ty,y(:,4),t);
weight = interp1(ty,y(:,5),t);
end

