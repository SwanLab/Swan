function strain = phitheta2strain(phi,theta)

stra_theta = theta;
stra_phi = phi;
stra_x = sin(stra_theta).*cos(stra_phi);
stra_y = sin(stra_theta).*sin(stra_phi);
stra_xy = cos(stra_theta);
strain = [stra_x stra_y stra_xy];

end