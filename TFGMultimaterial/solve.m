    %**************************************************************************
% Solve the pdetool FEM linear system  
%**************************************************************************
% DESCRIPTION
% Assemble and solve the pdetool fineite element linear system 
% 
% INPUT
% psi:     level-set function
% mesh:    pdetool mesh struct
% matprop: material properties struct
% pdecoef: pdetool coeficiets struct
% bc:      extra boundry conditions struct
%
% OUTPUT
% [U,F,vol]: pdetool solution, pdetool load and current volume
%
% HISTORY
% J-M.C. Farias     12/2010: code implementation.
% A.A. Novotny      12/2010: code updating.
%**************************************************************************

function [U,F,vol] = solve(psi, mesh, matprop, pdecoef, bc)
    
    gamma  = matprop.E./matprop.E(1); % contrast for each material
    p = mesh.p; e = mesh.e; t = mesh.t; area = mesh.area;
    b = pdecoef.b; a = pdecoef.a; f = pdecoef.f; c0 =  pdecoef.c(:,1);
    % element characteristic function: 
    [ ~,tfi ] = charfunc( p,t,psi);
%     tgamma = pdeintrp(p,t,fi*gamma');  %for P1-projection aproach
    tgamma = gamma*tfi; %for mixed formulation aproach
    c = c0 * tgamma; % effective elasticity tensor 
%     [K,F] = assempde(b,p,e,t,c,a,f);
    [K,~,F] = assema(p,t,c,a,f);
    [K,F] = pdeupdate(K,F,bc,mesh);
    U = K \ F;  % solve linear system
%     vol  = area*tchi';  %calculate volume assigned to each material including void -- P1 projection aproach
    vol  = area*tfi';  %calculate volume assigned to each material including void -- mixed formulation aproach
end

