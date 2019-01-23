function vademecum_mesh
close all
thetamin = 0;
thetamax = pi/2;
phimin = -pi/4;
phimax = pi/4;
cstring = 'rgbcmyk';

iteration = 1;
for ilevel = 1:15
    thetalevel = linspace(thetamin,thetamax,ilevel+1);
    
    for itheta = 1:length(thetalevel)
       philevel = linspace(phimin,phimax,itheta);
       if thetalevel(itheta) == 0
           philevel = 0;
       end
       
       for iphi = 1:length(philevel)
           phi = philevel(iphi);
           theta = thetalevel(itheta);
           deformation = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];          
           
           
           if ilevel == 1
           phi_all(iteration,1) = phi;
           theta_all(iteration,1) = theta;
           deformation_all(iteration,:) = deformation;
           
           figure(1)
           hold on
           plot(phi_all(iteration,1)*180/pi,theta_all(iteration,1)*180/pi,['+',cstring(mod(ilevel,7)+1)])
           hold off
           iteration = iteration + 1;
           end
           
           
           
           if ~((phi_all == phi) & (theta_all == theta) & ilevel > 1)
          
           phi_all(iteration,1) = phi;
           theta_all(iteration,1) = theta;
           deformation_all(iteration,:) = deformation;
  
           figure(1)
           hold on
           plot(phi_all(iteration,1)*180/pi,theta_all(iteration,1)*180/pi,['+',cstring(mod(ilevel,7)+1)])
           hold off
           iteration = iteration + 1;
           end
           
       end
    end
    figure(19)
    hold on
    plot3(deformation_all(:,1),deformation_all(:,2),deformation_all(:,3),['+',cstring(mod(ilevel,7)+1)])
    hold off
end    
end