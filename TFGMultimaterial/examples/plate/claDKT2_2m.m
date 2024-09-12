function [mesh, pdecoef, matprop, params, bc, psi] = claDKT2_2m
    %% geometry
    p = [ 0 , 0;
          10 , 0;
          10 , 10;
          0 , 10];
    g = [2; size(p,1); p(:,1); p(:,2)]; 

    % composing the geometry
    g = decsg(g);      
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
    params.penalty = [0.21];   % method 1,2 and 3
    params.volfrac = [0.6];   % method 2 and 3
    params.voleps  = [0.01];  % method 2 and 3
    params.auglag  = [5.1];   % method 3
    params.epsilon = 0.01;  % method 3      
    
    %% mesh generation
    n = 20;
    [p,e,t] = poimesh(g,n,n); t(4,:) = 1;
    ghold  = [];
    nsteps = 1;
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
    mesh.area = area;
    mesh.ghold = ghold;

    %% concentrated load-cases
    bc.pCons = []; %---> for static condensation
    
    bc.pNeu = [];  
    % lc.pNeu = [node,dof,val]
    % node: node number
    % dof:  degree of freedom number
    % val:  value of the nodal load 

    cx = 0.0; cy = 0.0;
    node = intersect(find(p(1,:)==cx),find(p(2,:)==cy));
    node = [node;node;node]; dof = [1;2;3]; val = [-1.0;0.0;0.0];
    bc.pNeu = [node, dof, val];
    
     
    %% Dirichlet boundary conditions
    
    bc.pDir = [];
    % bc.pDir = [node,dof,val]
    % node: node number
    % dof:  degree of freedom number
    % val:  value of the nodal constraint 

    
    cy = 10.0;
    node = find(p(2,:) == cy); nnodes = size(node,2); node = node';
    val0 = zeros(nnodes,1); dof1 = ones(nnodes,1); dof2 = 2*dof1; dof3 = 3*dof1;
    node = [node;node;node]; dof = [dof1;dof2;dof3]; val = [val0;val0;val0];
    pDir1 = [node, dof, val];
     
    cx = 10.0;
    node = find(p(1,:) == cx); nnodes = size(node,2); node = node';
    val0 = zeros(nnodes,1); dof1 = ones(nnodes,1); dof2 = 2*dof1; dof3 = 3*dof1;
    node = [node;node;node]; dof = [dof1;dof2;dof3]; val = [val0;val0;val0];
    pDir2 = [node, dof, val];
    
    cy = 0.0;
    node = find(p(2,:) == cy); nnodes = size(node,2); node = node';
    val0 = zeros(nnodes,1); dof2 = 2*ones(nnodes,1);
    node = [node]; dof = [dof2]; val = [val0];
    pDir3 = [node, dof, val];
     
    cx = 0.0;
    node = find(p(1,:) == cx); nnodes = size(node,2); node = node';
    val0 = zeros(nnodes,1); dof3 = 3*ones(nnodes,1);
    node = [node]; dof = [dof3]; val = [val0];
    pDir4 = [node, dof, val];
      
    bc.pDir = [pDir1; pDir2; pDir3; pDir4];

 
    %% Distributed Neumann boundary conditions    
    b = [];    
    
    %% material properties
    E0 = 210.0E+3;
    E = E0*[1 1e-4];
    h = 0.1;
    nu = 0.3;
    matprop.E = E;
    matprop.E0 = E;
    matprop.h0 = h;
    matprop.nu = nu; 
    matprop.nu0 = nu; 


        %% pde coeficients
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% FOR DKT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    c=zeros(2,length(matprop.E));
    c(1,:) = E*h*h*h./(12.*(1.-nu.*nu)); %Isotropic bending constitutive property
    c(2,:) = c(1)*nu;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% FOR DKMT & ZHUANG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%     k = 5/6;
%     c(3,:) = k*E*h./(2.*(1.+nu)); %Isotropic shear constitutive property
%     
%     
%     
%     pdecoef.c = c;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% FOR DKT & DKMT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pdecoef.a = zeros(9,1);
    pdecoef.f = zeros(3,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% FOR ZHUANG RM PLATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     pdecoef.a = zeros(15,1);
%     pdecoef.f = zeros(5,1);
    
    pdecoef.c = c;
    pdecoef.b = b;
    
        %% initial guess
    psi = ones(size(mesh.p,2),length(matprop.E)-1); psi(:,1) = -1;
end