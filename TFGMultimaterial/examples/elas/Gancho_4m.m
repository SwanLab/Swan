function [mesh, pdecoef, matprop, params, bc, psi] = viga2x1_4m

    %% load g and b
    load viga2x1.mat;
    
    %% parameters
  
    % minimum allowed 'k'. Used by the line-search procedure
    params.kmin = eps;
    % stop criterion
    params.stop = 1.0*pi/180;
    
    % method of penalization
    % (1) linear 
    % (2) exact quadratic
    % (3) augmented-lagrangian
    params.penalization = 3;
  
    % penalty parameter
    params.penalty = [0.0 0.0 0.0];   % method 1,2 and 3
    params.volfrac = [0.2 0.10 0.05];   % method 2 and 3
    params.voleps  = params.volfrac*0.05;  % method 2 and 3
    params.auglag  = [0.03 0.03 0.03];   % method 3
    params.epsilon = 0.01;  % method 3       
    
    
        
    %% mesh generation
    prec = 1e6;
    n = 10;
    [p,e,t] = poimesh(g,2*n,n); %p1 = round(prec*p1)/prec; 
    t(4,:) = 1; 
    
%     cd ..\..
%         [b,p,e,t] = raccommode(b1,p1,e1,t1,b2,p2,e2,t2);
%         [b,p,e,t] = raccommode(b,p,e,t,b3,p3,e3,t3);
% %         g = [g1,g2,g3];
%     cd examples\elas

    ghold  = 0;
    nsteps = 2;
    remesh  = 'longest';

    for i=1:2*nsteps
        [p,e,t] = refinemesh(g,p,e,t,remesh);            
    end

    area = pdetrg(p,t);

    mesh.remesh = remesh;
    mesh.p = p;
    mesh.e = e;
    mesh.t = t;
    mesh.g = g;
    mesh.area  = area;
    mesh.ghold = ghold;    
    
  %% extra boundary conditions
    
    bc.pCons = []; %---> for static condensation
  
    bc.pNeu = [];  
    % bc.pNeu = [node,dof,val]
    % node: node number
    % dof:  degree of freedom number
    % val:  value of the nodal load 
%-------------------------------------------------------------------------%
    bc.Neu.ldir = [];
    bc.Neu.ldof = [];
    bc.Neu.lval = [];

    ldir = [2.0,0.499,2.0,0.501];
        
    ldof = [1,2];
    
    lval = [-1];

    cd ..\..
    pbc = bouncon(p,e,t,ldir,ldof,lval);
    cd examples\elas
 
    bc.pNeu = pbc;        

    bc.Neu.ldir = ldir;
    bc.Neu.ldof = ldof;
    bc.Neu.lval = lval;  


    bc.pDir = [];
%-------------------------------------------------------------------------%    
%     % bc.pDir = [node,dof,val]
%     % node: node number
%     % dof:  degree of freedom number
%     % val:  value of the nodal constraint 
%-------------------------------------------------------------------------%
    ldir = [0.0,0.0,0.0,1.0];
    ldof = [2,1,2];
    lval = [0.0,0.0];
    
    cd ..\..
    
    pbc = bouncon(p,e,t,ldir,ldof,lval);      
      
    cd examples\elas
    
    bc.pDir = pbc;
    
    bc.Dir.ldir = ldir;
    bc.Dir.ldof = ldof;
    bc.Dir.lval = lval;
        %% Distributed Neumann boundary conditions    
    b = [];    
    %% material properties

    E = 200E9*[1 1/2 1/4 1/1000];     
    nu = 0.25;
    matprop.E = E;
    matprop.nu = nu; 

    %% pde coeficients
        
    la = nu.*E./((1+nu).*(1-2.*nu)); mu = E./(2.*(1+nu)); % plane strain
    la = 2.*mu.*la./(la+2.*mu); % plane stress
    
    c=zeros(16,length(E));
    
    c(1,:) = la + 2*mu; c(2,:) = 0; c(3,:) = 0; c(4,:) = mu;
    c(5,:) = 0; c(6,:) = la; c(7,:) = mu; c(8,:) = 0;
    c(9,:) = c(8,:); c(10,:) = c(7,:); c(11,:) = c(6,:); c(12,:) = c(5,:);
    c(13,:)= c(4,:); c(14,:) = c(3,:); c(15,:) = c(2,:); c(16,:) = c(1,:); 
    
    pdecoef.la = la;
    pdecoef.mu = mu;
    pdecoef.c = c;
    pdecoef.a = zeros(4,1);
    pdecoef.f = zeros(2,1);
    pdecoef.b = b;    
    
    %% initial guess
    psi = ones(size(mesh.p,2),length(matprop.E)-1); psi(:,1) = -1;
%     xc = 0.5; yc = 1; r = 0.25;
%     node =  (p(1,:)-xc).^2 + (p(2,:)-yc).^2 - r^2 <=0;
%     psi2(node) = -1;
%     psi1 = -rand(size(mesh.p,2),1)+0.5;
%     psi2 = -rand(size(mesh.p,2),1)+0.5;

end
