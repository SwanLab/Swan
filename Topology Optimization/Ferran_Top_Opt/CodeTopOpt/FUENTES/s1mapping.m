function [phi,theta] = s1mapping(eps_new)



   x = eps_new(1,:);
   y = eps_new(2,:);
   z = eps_new(3,:);
   
   for i = 1:length(eps_new(1,:))
   % Calculo de phi
   if x(i)> 0 && y(i)>=0
      phi(i,1) = atan(y(i)/x(i));
   end
    
   
   if x(i)> 0 && y(i)<0
      phi(i,1) = 2*pi + atan(y(i)/x(i));
   end
   
   
   
   if x(i)== 0 && y(i) < 0
      phi(i,1) = 3*pi/2;
   end
   
   if x(i)== 0 && y(i) >= 0
      phi(i,1) = pi/2;
   end
   
   
   
   if x(i)< 0
      phi(i,1) = pi + atan(y(i)/x(i));
   end
   
   
   
   if z(i)> 0
       theta(i,1) = atan(sqrt(x(i)^2+y(i)^2)/z(i));
   elseif z(i)==0
       theta(i,1) = pi/2;
   elseif z(i)<0
       theta(i,1) =  -atan(sqrt(x(i)^2+y(i)^2)/abs(z(i)));
   end
   
   
%     if phi(i,1) >= pi
%     phi(i,1) = phi(i,1) - pi;
%     end
    
    if theta(i,1) >= pi/2
    theta(i,1) = pi - theta(i,1);
    end

    phi = mod(phi,2*pi);
   end
end
   