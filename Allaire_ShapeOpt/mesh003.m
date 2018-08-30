function phi00 = mesh003(phi,RIiter,dx,dy,dz)

% REINITIALIZE THE LEVEL SET FUNCTION
% After intialization as well as periodically while
% solving the transport equation, we reinitialize
% the level set function, which means we solve 
% d(phi)/dt + sign(phi0)(|grad(phi)|-1) = 0 with
% phi(t=0,x) = phi_0(x)
% We do this to make sure the level set function
% does not get too flat or steep at the zero 
% level set, which is where our structure boundary
% is.  phi is the level set funtion we are
% reinitializing and RIiter is the number of times
% we're solving the reinitialization equation.       

  cfl = 0.5 ;   % CFL CONDITION
  % Using our CFL condition, we define the time step as
  dt0 = min(min(dx,dy),dz)*cfl ;

  for n = 1 : RIiter
    
    % For our first derivatives:
    phin = shift3n('n',phi,'n') ;
    phis = shift3n('s',phi,'n') ;
    phie = shift3n('e',phi,'n') ;
    phiw = shift3n('w',phi,'n') ;
    phiu = shift3n('u',phi,'n') ;
    phid = shift3n('d',phi,'n') ;

    % Our scheme is second order in space, so
    % we use these for our second derivatives:
    phinn = shift3n('n',phin,'n') ;
    phiss = shift3n('s',phis,'n') ;
    phiee = shift3n('e',phie,'n') ;
    phiww = shift3n('w',phiw,'n') ;
    phiuu = shift3n('u',phiu,'n') ;
    phidd = shift3n('d',phid,'n') ;
    
    
    %  The first derivatives:
    dxm = (phi-phie)/dx ;
    dxp = (phiw-phi)/dx ;
    dym = (phi-phis)/dy ;
    dyp = (phin-phi)/dy ;
    dzm = (phi-phid)/dz ;
    dzp = (phiu-phi)/dz ;
    
    
    
    % The second derivatives (our scheme is second
    % order so we use several different ones):
    dxmxm = (phi - 2*phie + phiee)/(dx^2) ;
    dxpxp = (phiww - 2*phiw + phi)/(dx^2) ;
    dxpxm = (phie - 2*phi + phiw)/(dx^2) ;

    dymym = (phi - 2*phis + phiss)/(dy^2) ;
    dypyp = (phinn - 2*phin + phi)/(dy^2) ;
    dypym = (phin - 2*phi + phis)/(dy^2) ;
    
    dzmzm = (phi - 2*phid + phidd)/(dz^2) ;
    dzpzp = (phiuu - 2*phiu + phi)/(dz^2) ;
    dzpzm = (phiu - 2*phi + phid)/(dz^2) ;
    
    
    
    % From Sethian, our scheme requires four
    % main parts, as defined here:
    partA = dxm + .5*dx*minmod(dxmxm,dxpxm) ;
    partB = dxp - .5*dx*minmod(dxpxp,dxpxm) ;
    partC = dym + .5*dy*minmod(dymym,dypym) ;
    partD = dyp - .5*dy*minmod(dypyp,dypym) ;
    
    partE = dzm + .5*dz*minmod(dzmzm,dzpzm) ;
    partF = dzp - .5*dz*minmod(dzpzp,dzpzm) ;
    
    
    
    
    % Then we use these parts along with the 
    % flux function to get:
    delp2 = g3(partA,partB,partC,partD,partE,partF) ;
    delm2 = g3(partB,partA,partD,partC,partF,partE) ;
    
    % sphi is the sign of phi and we added a small 
    % value proportional to the mesh size in the 
    % denominator to be sure we never divide by zero.
    nabla = 0.5*(dxm.^2 + dxp.^2 + dym.^2 + dyp.^2+ dzm.^2 + dzp.^2) ;
    sphi = phi./sqrt(phi.*phi+sqrt(dx^2+dy^2+dz^2)*nabla/10) ;
    sphip = max(sphi,0) ;
    sphim = min(sphi,0) ;

    phi = phi - dt0*(sphip.*delp2 + sphim.*delm2 - sphi) ;

  end

  % Once our iterations are finished, we set 
  % our output to the new level set function:
  phi00 = phi ;

end