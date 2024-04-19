function [mesh, pdecoef, matprop, params, bc, psi] = viga2x1_2m_40volfrac

    %% load g and b
    load viga2x1.mat;
    
    %% parameters
  
    % minimum allowed 'k'. Used by the line-search procedure
    params.kmin = eps;
    % stop criterion
    params.stop = 1.5*pi/180;
    
    % method of penalization
    % (1) linear 
    % (2) exact quadratic
    % (3) augmented-lagrangian
    params.penalization = 3;
  
    % penalty parameter
    params.penalty = [0.0];   % method 1,2 and 3
    params.volfrac = [0.4];   % method 2 and 3
    params.voleps  = [0.01];  % method 2 and 3
    params.auglag  = [1.7];   % method 3
    params.epsilon = 0.01;  % method 3       
    
    params.beta = 3.5; %perimetric constraint penalty parameter
    params.e_h = 32;  %perimetric constraint regularizing parameter ratio     
    %% mesh generation
    prec = 1e6;
    n = 80;
    [p,e,t] = poimesh(g,2*n,n); %p1 = round(prec*p1)/prec; 
    t(4,:) = 1; 
    
    params.scale.L = max(abs([min(min(p)) max(max(p))]));

    ghold  = 0;
    nsteps = 1;
    remesh  = 'longest';

    for i=1:2*nsteps
        [p,e,t] = refinemesh(g,p,e,t,remesh);            
    end
    p = p/params.scale.L;
    
    area = pdetrg(p,t);
    
    it1 = t(1,:); it2 = t(2,:); it3 = t(3,:);
    l12 = ((p(1,it1)-p(1,it2)).^2 + (p(2,it1)-p(2,it2)).^2).^(0.5);
    l13 = ((p(1,it1)-p(1,it3)).^2 + (p(2,it1)-p(2,it3)).^2).^(0.5);
    l32 = ((p(1,it3)-p(1,it2)).^2 + (p(2,it3)-p(2,it2)).^2).^(0.5);
    h = mean([l12 l13 l32]);

    mesh.remesh = remesh;
    mesh.p = p;
    mesh.e = e;
    mesh.t = t;
    mesh.g = g;
    mesh.area  = area;
    mesh.ghold = ghold;    
    mesh.h = h;
    
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
    
    xmax = max(p(1,:)); ymax = max(p(2,:));
    xmin = min(p(1,:)); ymin = min(p(2,:));
    
    ldir = [xmax,(ymax+ymin)/2 - eps,xmax,(ymax+ymin)/2 + eps];
        
    ldof = [1,2];
    
    lval = [-1];
    
    params.scale.t = max(abs(lval));

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
    
    %% material properties

    E = 70E9*[1 1/1000];   
    nu = 0.33;

    matprop.E = E;
    matprop.nu = nu; 
    matprop.gamma = E/max(E); 
    params.nmat = length(E);
    %% pde coeficients
        
    la = nu.*E./((1+nu).*(1-2.*nu)); mu = E./(2.*(1+nu)); % plane strain
    la = 2.*mu.*la./(la+2.*mu); % plane stress
    
    params.scale.mu = max(mu);

    mu = mu/params.scale.mu;
    la = la/params.scale.mu;
    
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
    %% Extra Scale Parameters
    params.scale.u_c = params.scale.t/params.scale.mu;
    
    %% initial guess
%     psi = ones(size(mesh.p,2),length(matprop.E)-1); psi(:,1) = -1;
    psi(:,1) = -(p(1,:)'-0.2).^2-(p(2,:)'-0.25).^2+0.1.^2;

end
