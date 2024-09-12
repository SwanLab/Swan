
%**************************************************************************
% Mesh Generator: the code have to be revised
%**************************************************************************
% DESCRIPTION
% Generates structured or symmetric meshes
% 
% INPUT
% g:        pdetool geometry
% nsteps:   number of mesh refinements
% meshtype: structured or regular mesh
% Hmax:     finite element size
%
% OUTPUT
% [p,e,t]:  pdetool mesh parameters
%
% HISTORY
% D.E. Campeão    12/2010: code implementation.
% A.A. Novotny  
%**************************************************************************

function [p,e,t] = genmesh(g,nsteps,meshtype,Hmax)
    % meshtype = 'longest' or 'regular'
    [p,e,t] = initmesh(g,'Hmax',Hmax,'Init','on');
    [p,e,t] = refinemesh(g,p,e,t,'longest');
    p = jigglemesh(p,e,t,'Opt','minimum','Iter',100);
    for i=1:(strcmp(meshtype,'longest')*2+strcmp(meshtype,'regular'))*nsteps
        [p,e,t] = refinemesh(g,p,e,t,meshtype);            
    end
    p = jigglemesh(p,e,t,'Opt','minimum','Iter',100);
    pdemesh(p,e,t);
end