function [mesh, pdecoef, matprop, params, bc, psi] = star_4m

    %% load g and b
    load star.mat;
    p = [ 0 , 0;
          1 , 0;
          1 , 1;
          0 , 1];
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
    params.penalty = [0.0 0.0 0.0];   % method 1,2 and 3
    params.volfrac = [0.25 0.25 0.25];   % method 2 and 3
    params.voleps  = [0.01 0.01 0.01];  % method 2 and 3
    params.auglag  = [0.02 0.02 0.02];   % method 3
    params.epsilon = 0.01;  % method 3     
        
    %% mesh generation
    n = 40;
    [p,e,t] = poimesh(g,n,n); 
    
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
    mesh.area  = area;
    mesh.ghold = ghold;    
    
    %% extra boundary conditions
    bc.pCons = []; %---> for static condensation
    
    bc.pNeu = [];  
    % bc.pNeu = [node,dof,val]
    % node: node number
    % dof:  degree of freedom number
    % val:  value of the nodal load 
    bc.Neu.ldir = [];
    bc.Neu.ldof = [];
    bc.Neu.lval = [];

    cx = [0.0,0.3;...
          0.7,1.0];
    cy = [1.0,1.0;...
          0.0,0.0];

    ldir = [cx(:,1),cy(:,1),cx(:,2),cy(:,2)];
    ldof = ones(2,2);
    lval = -ones(2,1);

    cd ..\..
        pbc = bouncon(p,e,t,ldir,ldof,lval);
    cd examples\elas
    
    bc.pNeu = pbc;
        
%-------------------------------------------------------------------------%    
%     % bc.pDir = [node,dof,val]
%     % node: node number
%     % dof:  degree of freedom number
%     % val:  value of the nodal constraint 
%-------------------------------------------------------------------------%
    cx = [0.0,0.0;...
          1.0,1.0]; 
    cy = [0.0,0.3;...
          0.7,1.0];

    ldir = [cx(:,1),cy(:,1),cx(:,2),cy(:,2)];
    ldof = ones(2,2);
    lval = zeros(2,1);
    
    pDir = [];
    cd ..\..
        pbc = bouncon(p,e,t,ldir,ldof,lval);
    cd examples\elas
    
    bc.pDir = pbc;
    
    bc.Dir.ldir = ldir;
    bc.Dir.ldof = ldof;
    bc.Dir.lval = lval;
        %% Distributed Neumann boundary conditions    
%     b = [];    
    %% material properties
    
    k = 1e3.*[1 0.5 0.25 1e-4];     
    matprop.k = k;
    matprop.E = k;
    
    %% pde coeficients
      
    c = k;
    pdecoef.c = c;
    pdecoef.a = 0.0;
    pdecoef.f = 0.0;
    pdecoef.b = b;    
    
    %% initial guess
    psi = ones(size(mesh.p,2),length(matprop.k)-1); psi(:,1) = -1;

end
