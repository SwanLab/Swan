function plot_gp_field(field,coordinates,dim,problembsc,element)
xgp = interpol(coordinates(:,1),element,dim,problembsc,coordinates);
ygp = interpol(coordinates(:,2),element,dim,problembsc,coordinates);

tri = delaunay(xgp,ygp);
figure
trisurf(tri,xgp,ygp,field);

end