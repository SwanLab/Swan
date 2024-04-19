function [mesh, pdecoef, matprop, params, bc, psi] = viga2x1

    %% load g and b
    load viga2x1.mat;
    
    %% parameters
  
    % minimum allowed 'k'. Used by the line-search procedure
    params.kmin = 1.0E-4;
    % stop criterion
    params.stop = 1.0*pi/180;
    
    % method of penalization
    % (1) linear 
    % (2) exact quadratic
    % (3) augmented-lagrangian
    params.penalization = 3;
  
    % penalty parameter
    params.penalty(1) = 0.0;   % method 1,2 and 3
    params.penalty(2) = 0.0;   % method 1,2 and 3
    params.volfrac(1) = 0.15;   % method 2 and 3
    params.volfrac(2) = 0.15;   % method 2 and 3
    params.voleps(1)  = 0.01;  % method 2 and 3
    params.voleps(2)  = 0.01;  % method 2 and 3
    params.auglag(1)  = 2;   % method 3
    params.auglag(2)  = 2;   % method 3
    params.epsilon = 0.01;  % method 3       
    
    
        
    %% mesh generation
    prec = 1e6;
    n = 8;
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
    
    %% material properties

    E1 =200E9; 
    E2 =50E9;
    E3 = E1*0.0001; 
    nu = 0.3;
    
    matprop.E(1) = E1;
    matprop.E(2) = E2; 
    matprop.E(3) = E3; 
    matprop.nu = nu; 

    %% pde coeficients
    
    %     Material 1
    la1 = nu*E1/((1+nu)*(1-2*nu)); mu1 = E1/(2*(1+nu)); % plane strain
    la1 = 2*mu1*la1/(la1+2*mu1); % plane stress
    
    c1=zeros(16,1);
    c1(1) = la1 + 2*mu1; c1(2) = 0; c1(3) = 0; c1(4) = mu1;
    c1(5) = 0; c1(6) = la1; c1(7) = mu1; c1(8) = 0; 
    c1(9) = c1(8); c1(10) = c1(7); c1(11) = c1(6); c1(12) = c1(5);
    c1(13)= c1(4); c1(14) = c1(3); c1(15) = c1(2); c1(16) = c1(1);

    %     Material 2
    la2 = nu*E2/((1+nu)*(1-2*nu)); mu2 = E2/(2*(1+nu)); % plane strain
    la2 = 2*mu2*la2/(la2+2*mu2); % plane stress
    
    c2=zeros(16,1);
    c2(1) = la2 + 2*mu2; c2(2) = 0; c2(3) = 0; c2(4) = mu2;
    c2(5) = 0; c2(6) = la2; c2(7) = mu2; c2(8) = 0; 
    c2(9) = c2(8); c2(10) = c2(7); c2(11) = c2(6); c2(12) = c2(5);
    c2(13)= c2(4); c2(14) = c2(3); c2(15) = c2(2); c2(16) = c2(1);
        
    %     Material 3  
    la3 = nu*E3/((1+nu)*(1-2*nu)); mu3 = E3/(2*(1+nu)); % plane strain
    la3 = 2*mu3*la3/(la3+2*mu3); % plane stress
  
    c3=zeros(16,1);
    c3(1) = la3 + 2*mu3; c3(2) = 0; c3(3) = 0; c3(4) = mu3;
    c3(5) = 0; c3(6) = la3; c3(7) = mu3; c3(8) = 0; 
    c3(9) = c3(8); c3(10) = c3(7); c3(11) = c3(6); c3(12) = c3(5);
    c3(13)= c3(4); c3(14) = c3(3); c3(15) = c3(2); c3(16) = c3(1);
    
    
    pdecoef.la(1) = la1;
    pdecoef.mu(1) = mu1;
    pdecoef.la(2) = la2;
    pdecoef.mu(2) = mu2;
    pdecoef.la(2) = la3;
    pdecoef.mu(2) = mu3;
    pdecoef.c(:,1) = c1;
    pdecoef.c(:,2) = c2;
    pdecoef.c(:,3) = c3;
    pdecoef.a = zeros(4,1);
    pdecoef.f = zeros(2,1);
    pdecoef.b = b;
    
    
    %% initial guess
    psi1 = -ones(size(mesh.p,2), 1);
    psi2 = ones(size(mesh.p,2), 1);
%     xc = 0.5; yc = 1; r = 0.25;
%     node =  (p(1,:)-xc).^2 + (p(2,:)-yc).^2 - r^2 <=0;
%     psi2(node) = -1;
%     psi1 = -rand(size(mesh.p,2),1)+0.5;
%     psi2 = -rand(size(mesh.p,2),1)+0.5;
    psi = [psi1 psi2];

end
