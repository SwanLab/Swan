%**************************************************************************
% solution update   
%**************************************************************************
% DESCRIPTION
% Updates solution of the linear system
% 
% INPUT
% U:    pdetool reduced solution
% bc:   extra boundary conditions struct
%
% OUTPUT
% [U,F]: updated pdetool solution and load vector
%
% HISTORY
% S.M. Giusti     02/2011: code implementation.
%
%**************************************************************************

function [U,F] = solupdate(Ur,bc,mesh,a,b,c,f)

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
        nbdofa = length(bc.pCons.dofs.all);
        nbdofi = length(bc.pCons.dofs.i);
        nbdofp = length(bc.pCons.dofs.p);
        nbdofm = length(bc.pCons.dofs.m);
        
        U = zeros(nbdofa,1);
        
        U(bc.pCons.dofs.i,:) = Ur(1:nbdofi,:);
        U(bc.pCons.dofs.p,:) = Ur(nbdofi+1:length(Ur),:);
        U(bc.pCons.dofs.m,:) = Ur(nbdofi+1:length(Ur),:);
        
        % Retrieve original forcing term
        p = mesh.p; e = mesh.e; t = mesh.t;
        [~,F] = assempde(b,p,e,t,c,a,f);
        aa=1
    else
        U = Ur;
        % Retrieve original forcing term
        p = mesh.p; e = mesh.e; t = mesh.t;
        [~,F] = assempde(b,p,e,t,c,a,f);
        aa=2
    end
    

end