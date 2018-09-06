function solvedphi = solvelvlset3(phi,V,dt,HJiter,lagP,RIiter,RIfreq,dx,dy,dz)

% SOLVE THE TRANSPORT LEVEL SET EQUATION
% This function solves an Hamilton Jacobi equation
% to transport the level set function:
% d(phi)/dt + V |grad(phi)| = 0 with
% where V is the normal velocity (the shape gradient
% of the objective function).
% The scheme is second order in space and explicit
% first order in time. The time step is dt and the
% number of time steps (or iterations) is HJiter.

% figure('NumberTitle', 'off', 'Name', 'FEM-MAT-OO - Phi_n')
% subplot(2,2,1), surf(phi(:,:,1)), title('\phi_n - Root')
% subplot(2,2,2), surf(phi(:,:,end)), title('\phi_n - Tip')
% subplot(2,2,3), surf(permute(phi(ceil(size(phi,1)/2),:,:),[2 3 1])), title('\phi_n - XY')
% subplot(2,2,4), surf(permute(phi(:,ceil(size(phi,2)/2),:),[1 3 2])), title('\phi_n - XZ')

for i = 1 : HJiter
    
    % We reinitialize the level set function for every RIfreq
    % time steps while solving the transport level set equation:
    if mod(i,RIfreq) == 0 || i == 1
        phi = mesh003(phi,RIiter,dx,dy,dz);
    end
    
    % Our scheme takes into account whether our
    % front is moving forward or backward:
    Vp = max(V,0) ;
    Vm = min(V,0) ;
    
    % We use these for our first derivatives:
    phin = shift3n('n',phi,'ng') ;
    phis = shift3n('s',phi,'ng') ;
    phie = shift3n('e',phi,'ng') ;
    phiw = shift3n('w',phi,'ng') ;
    phiu = shift3n('u',phi,'ng') ;
    phid = shift3n('d',phi,'ng') ;
    
    
    % We use these for our second derivatives:
    
    phinn = shift3n('n',phin,'ng') ;
    phiss = shift3n('s',phis,'ng') ;
    phiee = shift3n('e',phie,'ng') ;
    phiww = shift3n('w',phiw,'ng') ;
    phiuu = shift3n('u',phiu,'ng') ;
    phidd = shift3n('d',phid,'ng') ;
    
    
    % The first derivatives:
    dxm = (phi-phie)/dx ;
    dxp = (phiw-phi)/dx ;
    dym = (phi-phis)/dy ;
    dyp = (phin-phi)/dy ;
    dzm = (phi-phid)/dz ;
    dzp = (phiu-phi)/dz ;
    
    % Because our scheme is second order in space,
    % we use the second derivatives:
    dxmxm = (phi - 2*phie + phiee)/(dx^2) ;
    dxpxp = (phiww - 2*phiw + phi)/(dx^2) ;
    dxpxm = (phie - 2*phi + phiw)/(dx^2) ;
    
    dymym = (phi - 2*phis + phiss)/(dy^2) ;
    dypyp = (phinn - 2*phin + phi)/(dy^2) ;
    dypym = (phin - 2*phi + phis)/(dy^2) ;
    
    dzmzm = (phi - 2*phid + phidd)/(dz^2) ;
    dzpzp = (phiuu - 2*phiu + phi)/(dz^2) ;
    dzpzm = (phiu - 2*phi + phid)/(dz^2) ;
    
    % Our scheme uses certain parts to which
    % we apply the flux function
    partA = dxm + .5*dx*minmod(dxmxm,dxpxm) ;
    partB = dxp - .5*dx*minmod(dxpxp,dxpxm) ;
    partC = dym + .5*dy*minmod(dymym,dypym) ;
    partD = dyp - .5*dy*minmod(dypyp,dypym) ;
    
    partE = dzm + .5*dz*minmod(dzmzm,dzpzm) ;
    partF = dzp - .5*dz*minmod(dzpzp,dzpzm) ;
    
    delp2 = g3(partA,partB,partC,partD,partE,partF) ;
    delm2 = g3(partB,partA,partD,partC,partF,partE) ;
    
    % Here we compute the |grad(phi)|.
    
    if lagP==0
        % We update our old level set:
        phi = phi - dt*(delp2.*Vp + delm2.*Vm);
        
    else
        epscurv =min( min(dx,dy),dz)/20 ;
        
        dphix = (phiw-phie)/(2*dx) ;
        dphiy = (phin-phis)/(2*dy) ;
        dphiz = (phiu-phid)/(2*dz) ;
        
        mag = sqrt(dphix.^2+dphiy.^2+dphiz.^2+epscurv^2) ;
        
        % We update our old level set:
        phi = phi - dt*(delp2.*Vp + delm2.*Vm)+dt*lagP*curv3(phi,dx,dy,dz).*mag ;
    end
    
end

% We let our solved level set be the output:
solvedphi = phi ;
%   phi(obj.FEM.row2mat( obj.unrem_nodes_index,'nodes')==1)=-0.1;

% figure, surf(-permute(phi(:,ceil(size(phi,2)/2),:),[1 3 2])), title('phi')

% figure('NumberTitle', 'off', 'Name', 'Jesus - Phi_n_+_1')
% subplot(2,2,1), surf(phi(:,:,1)), title('\phi_n_+_1 - Root')
% subplot(2,2,2), surf(phi(:,:,end)), title('\phi_n_+_1 - Tip')
% subplot(2,2,3), surf(permute(phi(ceil(size(phi,1)/2),:,:),[2 3 1])), title('\phi_n_+_1 - XY')
% subplot(2,2,4), surf(permute(phi(:,ceil(size(phi,2)/2),:),[1 3 2])), title('\phi_n_+_1 - XZ')
end
