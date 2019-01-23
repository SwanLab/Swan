[Coordinates,Elements] = callesferita;
ncoord = size(Coordinates,1);
[phi,theta] = s1mapping([Coordinates(:,2) Coordinates(:,3) zeros(ncoord,1) Coordinates(:,4)]');

phi_logical = phi >= (360 -45)*pi/180;
phi(phi_logical) = phi(phi_logical) - 2*pi;

nodes = (phi <= pi/4 & phi >= -pi/4 & theta <= pi/2);

plot(phi(nodes)*180/pi,theta(nodes)*180/pi,'+')