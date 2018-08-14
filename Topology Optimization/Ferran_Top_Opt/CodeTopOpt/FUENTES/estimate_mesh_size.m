function  h = estimate_mesh_size(element,coordinates)
x1 = coordinates(element.conectivities(:,1));
x2 = coordinates(element.conectivities(:,2));
x3 = coordinates(element.conectivities(:,3));

x1x2 = abs(x2 - x1);
x2x3 = abs(x3 - x2);
x1x3 = abs(x1-x3);
hs = max([x1x2,x2x3,x1x3]');
h = mean(hs);


end