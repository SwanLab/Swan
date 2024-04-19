%**************************************************************************
% pdetool update   
%**************************************************************************
% DESCRIPTION
% Updates pdetool linear system
% 
% INPUT
% K:    pdetool stiffness matrix
% F:    pdetool load vector
% bc:   extra boundary conditions struct
% mesh: mesh parameters struct 
%
% OUTPUT
% [K,F]: updated pdetool stiffness matrix and load vector
%
% HISTORY
% J-M.C. Farias     12/2010: code implementation.
% D.E. Campeão     
% A.A. Novotny
%**************************************************************************

function [K,F] = pdeupdate(K,F,bc,mesh)

    np = size(mesh.p,2);   
    if ~isempty(bc.pNeu)
       %bc.pNeu=[node,dof,value]
        eq = bc.pNeu(:,1) + (bc.pNeu(:,2)-1)*np;
        F(eq) = F(eq) + bc.pNeu(:,3); 
    end

    if ~isempty(bc.pDir)
        %bc.pDir=[node,dof,value]
        eq = bc.pDir(:,1) + (bc.pDir(:,2)-1)*np ; 
        U = zeros(size(F,1),1);
        U(eq) = bc.pDir(:,3);
        F = F - K * U;
        K(eq,:)=0.0 ;    K(:,eq)=0.0 ;
        for i=eq'
            K(i,i)=1.0;
        end
        F(eq) = bc.pDir(:,3);
    end
    
    if ~isempty(bc.pCons)
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
    
        % Set relevant stiffness sub-matrices
        Kii = K(bc.pCons.dofs.i(:),bc.pCons.dofs.i(:));
        Kip = K(bc.pCons.dofs.i(:),bc.pCons.dofs.p(:));
        Kim = K(bc.pCons.dofs.i(:),bc.pCons.dofs.m(:));
        Kpi = K(bc.pCons.dofs.p(:),bc.pCons.dofs.i(:));
        Kpp = K(bc.pCons.dofs.p(:),bc.pCons.dofs.p(:));
        Kpm = K(bc.pCons.dofs.p(:),bc.pCons.dofs.m(:));
        Kmi = K(bc.pCons.dofs.m(:),bc.pCons.dofs.i(:));
        Kmp = K(bc.pCons.dofs.m(:),bc.pCons.dofs.p(:));
        Kmm = K(bc.pCons.dofs.m(:),bc.pCons.dofs.m(:));
        
        % Assemble reduced system matrix
        Kr = [   Kii         Kip+Kim      ;
               Kpi+Kmi   Kpp+Kpm+Kmp+Kmm ];
        
        % Set relevant forcing term sub-vectors
        Fi = F(bc.pCons.dofs.i(:));
        Fp = F(bc.pCons.dofs.p(:));
        Fm = F(bc.pCons.dofs.m(:));
        
        % Assemble reduced system forcing term
        Fr = [  Fi  ;
              Fp+Fm];
        
        K = Kr;
        F = Fr;

    end  
    
end