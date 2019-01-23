% THE PERIMETER
% In order to roughly calculate the perimeter
% we find the norm of the gradient of the sign
% of phi and integrate it, then divide by 2:
function totperim = perimeter(phi,dx,dy)
  % To smooth sign(phi):
  epsperim = min(dx,dy)/20 ;
  sx = phi./sqrt(phi.^2+epsperim^2) ;
    
  sxn = shift2n('n',sx,'ng') ;
  sxs = shift2n('s',sx,'ng') ;
  sxe = shift2n('e',sx,'ng') ;
  sxw = shift2n('w',sx,'ng') ;
  
  % We now calculate d(phi)/dx and d(phi)/dy:
  dsxx = (sxw-sxe)/(2*dx) ;
  dsxy = (sxn-sxs)/(2*dy) ;
  
  dV = dx*dy ;
  
  % And then integrate:
  totperim = .5*sum(sum(sqrt(dsxx.^2+dsxy.^2)))*dV ;
end
