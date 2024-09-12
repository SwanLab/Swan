function [mesh, pdecoef, matprop, params, bc, psi] = cog_4m
    %% load g and b 
        %% Distributed Neumann boundary conditions    
    load cog.mat
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
    params.penalty = [0.0 0.0 0.0];   % method 1,2 and 3
    params.volfrac = 0.5*[1/3 1/3 1/3];   % method 2 and 3
    params.voleps  = [0.01 0.01 0.01];  % method 2 and 3
    params.auglag  = [0.02 0.02 0.02];   % method 3
    params.epsilon = 0.01;  % method 3      
   
    %% mesh generation  
    n = 40;
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
    
    cx = 0.0; cy = 0.0;    
    ldir = [cx,cy-2*eps,cx,cy+2*eps];        
    ldof = [1,1];    
    lval = [-1];
    
    cd ..\..
    pbc = bouncon(p,e,t,ldir,ldof,lval);
    cd examples\elas
 
    bc.pNeu = pbc;        

    bc.Neu.ldir = ldir;
    bc.Neu.ldof = ldof;
    bc.Neu.lval = lval;  

    bc.pDir = [];
    % bc.pDir = [node,dof,val]
    % node: node number
    % dof:  degree of freedom number
    % val:  value of the nodal constraint 

    cx = max(p(1,:)); cy = max(p(2,:));
    ldir1 = [0,cy,cx,cy];
    ldir2 = [cx,0,cx,cy];
    ldir = [ldir1;ldir2];
    ldof = [3,1,2,3;...
            3,1,2,3];
    lval = [0.0,0.0,0.0;...
            0.0,0.0,0.0];
    
    cd ..\..
    
    pbc = bouncon(p,e,t,ldir,ldof,lval);      
      
    cd examples\elas
    
    bc.pDir = pbc;
    
    bc.Dir.ldir = ldir;
    bc.Dir.ldof = ldof;
    bc.Dir.lval = lval;
    

    %% material properties
    E0 = 210.0E+3;
    E = E0*[1 0.5 0.25 1e-4];
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
    
    
    
    pdecoef.c = c;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% FOR DKT & DKMT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pdecoef.a = zeros(9,1);
    pdecoef.f = zeros(3,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% FOR ZHUANG RM PLATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     pdecoef.a = zeros(15,1);
%     pdecoef.f = zeros(5,1);
    
    pdecoef.c = c;
    pdecoef.b = b;
    
    
    %% initial guess
    psi = ones(size(mesh.p,2),length(matprop.E)-1); psi(:,1) = -1;
    
end

