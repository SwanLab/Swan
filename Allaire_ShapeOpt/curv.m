
% THE CURVATURE
% The shape gradient of the perimeter is the  
% mean curvature which we compute by 
% div(grad(phi)/|grad(phi)|)
% where phi is the level set function. 
function H = curv(phi,dx,dy)
  % When finding the normal vector, we need
  % to divide by the norm of the gradient of phi;
  % we use this small value to make sure that
  % the gradient of phi never goes to 0:
  epscurv = min(dx,dy)/20 ;
  
  % Here are the first derivatives for finding
  % the gradient of phi:
  phin = shift2n('n',phi,'ng') ;
  phis = shift2n('s',phi,'ng') ;
  phie = shift2n('e',phi,'ng') ;
  phiw = shift2n('w',phi,'ng') ;

  dphix = (phiw-phie)/(2*dx) ;
  dphiy = (phin-phis)/(2*dy) ;
  
  %  Here is |grad(phi)| and then the x and y
  % components of the normal vector field:
  mag = sqrt(dphix.^2+dphiy.^2+epscurv^2) ;
  nx = dphix./mag ; ny = dphiy./mag ;
  
  % Now to find the divergence, just take the 
  % partials with repsect to x and y of the
  % x and y components of our normal vector field
  % (respectively) and then add these together
  % to get our end mean curvature, which is a 
  % function across our working domain.
  nxe = shift2n('e',nx,'ng') ;
  nxw = shift2n('w',nx,'ng') ;
  
  nyn = shift2n('n',ny,'ng') ;
  nys = shift2n('s',ny,'ng') ;
  
  divnx = (nxw-nxe)/(2*dx) ;
  divny = (nyn-nys)/(2*dy) ;
    
  H = divnx+divny ;
end