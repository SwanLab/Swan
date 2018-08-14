function [shape] = shape_func_triang(dirichlet_data,coord,e,x,y)
% e : sublista de elementos
% (x,y) : puntos 

shape = zeros(3,size(e,1));
x1 = coord(dirichlet_data(1,e),1);
x2 = coord(dirichlet_data(2,e),1);
x3 = coord(dirichlet_data(3,e),1);
y1 = coord(dirichlet_data(1,e),2);
y2 = coord(dirichlet_data(2,e),2);
y3 = coord(dirichlet_data(3,e),2);
area2 = (x2.*y3-x3.*y2) - (x1.*y3-x3.*y1) + (x1.*y2-x2.*y1);
shape(1,:)=(x2.*y3 - x3.*y2 + (y2-y3).*x + (x3-x2).*y)./area2;
shape(2,:)=(x3.*y1 - x1.*y3 + (y3-y1).*x + (x1-x3).*y)./area2;
shape(3,:)=(x1.*y2 - x2.*y1 + (y1-y2).*x + (x2-x1).*y)./area2;
end

