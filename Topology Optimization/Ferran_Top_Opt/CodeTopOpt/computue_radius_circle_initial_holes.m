%Compute radius circle
Vfrac = 0.6;
ncircles = 3;
geometric_domain = 2;
circle_area = (1 - Vfrac)*geometric_domain/ncircles;
r = sqrt(circle_area/pi)
length = 2;
h = length/ncircles;
center_position = h/2:h:2