function plot_nodal_field(field,coordinates)
x = coordinates(:,1);
y = coordinates(:,2);

tri = delaunay(x,y);
trisurf(tri,x,y,field);

end