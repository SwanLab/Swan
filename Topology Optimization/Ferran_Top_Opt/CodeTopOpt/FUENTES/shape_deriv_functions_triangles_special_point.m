function  shape = shape_deriv_functions_triangles_special_point(posgp)
s=posgp(1,:);t=posgp(2,:); 
%Triangle
shape(1,:)=1.0-s-t;
shape(2,:)=s;
shape(3,:)=t;



end