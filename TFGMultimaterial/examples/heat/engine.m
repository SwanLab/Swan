function [mesh, pdecoef, matprop, params, bc, psi0] = engine

    %% load g and b
    % load engine.mat;
     load engine2.mat;
    % load engine3.mat;
    
    %% parameters
  
    % minimum allowed 'k'. Used by the line-search procedure
    params.kmin = 1.0E-3;
    % stop criterion
    params.stop = 1.0*pi/180;
    
    % method of penalization
    % (1) linear 
    % (2) exact quadratic
    % (3) augmented-lagrangian
    params.penalization = 3;
  
    % penalty parameter
    params.penalty = 0.6;  % method 1, 2 and 3
    params.volfrac = 0.635;  % method 2 and 3
    params.voleps  = 0.01; % method 2 and 3
    params.auglag  = 0.1;  % method 3
    params.epsilon = 0.01; % method 3       
        
    %% mesh generation
    [p,e,t] = poimesh(g,20,20); t(4,:) = 1;
    
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
    
    bc.pNeu = [];  
    % bc.pNeu = [node,dof,val]
    % node: node number
    % dof:  degree of freedom number
    % val:  value of the nodal load 

    bc.pDir = [];
    % bc.pDir = [node,dof,val]
    % node: node number
    % dof:  degree of freedom number
    % val:  value of the nodal constraint 
        

    bc.pCons = [];
    % bc.pCon = bc.pCons.node
    %           bc.pCons.node.i
    %           bc.pCons.node.p
    %           bc.pCons.node.m
    %           bc.pCons.node.all
    % bc.pCon = bc.pCons.dofs
    %           bc.pCons.dofs.i
    %           bc.pCons.dofs.p
    %           bc.pCons.dofs.m
    %           bc.pCons.dofs.all

    % node: node number
    % dof:  degree of freedom number
     bc.pCons.node = [];
     bc.pCons.dofs = [];

    

    %% material properties
    
    k = 1.0; k1 = 1.0E-5;
    matprop.k0 = k; matprop.k1 = k1;
    matprop.E0 = k; matprop.E1 = k1;
    
    %% pde coeficients
      
    c = k;
    pdecoef.c = c;
    pdecoef.a = 0.0;
    pdecoef.f = 0.0;
    pdecoef.b = b;    
    
    %% initial guess
    psi0 = -ones(size(mesh.p,2), 1);

end
