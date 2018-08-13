
function solvedphi = solvelvlset(phi,V,dt,HJiter,lagP,RIiter,RIfreq,dx,dy)
  
% SOLVE THE TRANSPORT LEVEL SET EQUATION 
% This function solves an Hamilton Jacobi equation 
% to transport the level set function:
% d(phi)/dt + V |grad(phi)| = 0 with
% where V is the normal velocity (the shape gradient 
% of the objective function).
% The scheme is second order in space and explicit 
% first order in time. The time step is dt and the 
% number of time steps (or iterations) is HJiter.
  for i = 1 : HJiter 
    
    % We reinitialize the level set function for every RIfreq
    % time steps while solving the transport level set equation:
    if mod(i,RIfreq) == 0 || i == 1
      phi = mesh00(phi,RIiter,dx,dy) ;
    end     
    
    % Our scheme takes into account whether our
    % front is moving forward or backward:
    Vp = max(V,0) ;
    Vm = min(V,0) ;
    
    % We use these for our first derivatives:
    phin = shift2n('n',phi,'ng') ;
    phis = shift2n('s',phi,'ng') ;
    phie = shift2n('e',phi,'ng') ;
    phiw = shift2n('w',phi,'ng') ;
    
    % We use these for our second derivatives:
    phinn = shift2n('n',phin,'ng') ;
    phiss = shift2n('s',phis,'ng') ;
    phiee = shift2n('e',phie,'ng') ;
    phiww = shift2n('w',phiw,'ng') ;
    
    % The first derivatives:
    dxm = (phi-phie)/dx ;
    dxp = (phiw-phi)/dx ;
    dym = (phi-phis)/dy ;
    dyp = (phin-phi)/dy ;
    
    % Because our scheme is second order in space,
    % we use the second derivatives:
    dxmxm = (phi - 2*phie + phiee)/(dx^2) ;
    dxpxp = (phiww - 2*phiw + phi)/(dx^2) ;
    dxpxm = (phie - 2*phi + phiw)/(dx^2) ;

    dymym = (phi - 2*phis + phiss)/(dy^2) ;
    dypyp = (phinn - 2*phin + phi)/(dy^2) ;
    dypym = (phin - 2*phi + phis)/(dy^2) ;
    
    % Our scheme uses certain parts to which
    % we apply the flux function
    partA = dxm + .5*dx*minmod(dxmxm,dxpxm) ;
    partB = dxp - .5*dx*minmod(dxpxp,dxpxm) ;
    partC = dym + .5*dy*minmod(dymym,dypym) ;
    partD = dyp - .5*dy*minmod(dypyp,dypym) ;

    delp2 = g(partA,partB,partC,partD) ;
    delm2 = g(partB,partA,partD,partC) ;
    
    % Here we compute the |grad(phi)|.
    epscurv = min(dx,dy)/20 ;

    dphix = (phiw-phie)/(2*dx) ;
    dphiy = (phin-phis)/(2*dy) ;
  
    mag = sqrt(dphix.^2+dphiy.^2+epscurv^2) ;
    
    % We update our old level set:
    phi = phi - dt*(delp2.*Vp + delm2.*Vm)+dt*lagP*curv(phi,dx,dy).*mag ;
    
%     surf(phi), view([0 0 1])
  end 
  
  % We let our solved level set be the output:
  solvedphi = phi ;
  
end

