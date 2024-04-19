%**************************************************************************
% Find nodes on the boundary of the MEF on an specific location
%**************************************************************************
%
% DESCRIPTION
% Find nodes on the boundary of the MEF to apply boundary conditions.
%
% HISTORY
% S.M. Giusti      07/2012: code implementation.
%
% INPUT
% p      : mesh data from pdetool (points)
% e      : mesh data from pdetool (edges)
% t      : mesh data from pdetool (triangules)
% lnodes : lines to find nodes. The lines are defined in a row with four
%          parameters: coordinates of the starting and ending points of
%          each line. Example:
%                         lnodes = [xi,yi,xf,yf]   or 
%                         lnodes = [[ 5.0,0.0,10.0,0.0];...
%                                   [10.0,0.0,10.0,3.0]]
% OUTPUT
% nodes : nodes in the lnodes array
% edges : edges of the FEM mesh in lnodes array
%**************************************************************************

function [nodes,edges,triangule] = findnodes(p,e,t,lnodes)
    
    nlines = size(lnodes,1);
    nedge = e(1:2,:);
    
    nodes = [];
    
    for i=1:nlines
        
        iline = lnodes(i,:);
        
        cxi = iline(1);  cyi = iline(2);
        cxf = iline(3);  cyf = iline(4);
        
        nodex = intersect(find(p(1,:) >= cxi),find(p(1,:) <= cxf));
        nodey = intersect(find(p(2,:) >= cyi),find(p(2,:) <= cyf));
        
        node = intersect(nodex,nodey);
        node = intersect(node,unique(e(1:2,:)))';
        nodes = [nodes;node];
        
    end
    nodes = unique(nodes);
    
    edgen1 = [];
    edgen2 = [];

    for i=1:size(nodes,1)
        
        aux = find(nedge(1,:) == nodes(i));
        if ~isempty(aux)
            %edgen1(i) = aux;
            edgen1 = [edgen1;aux'];
        end
        aux = find(nedge(2,:) == nodes(i));
        if ~isempty(aux)
            %edgen2(i) = aux;
            edgen2 = [edgen2;aux'];
        end
    end
    
    edges=intersect(edgen1,edgen2);
    
    if isempty(edges)
        triangule = [];
        edges = [];
    else
        conedges = nedge(:,edges);
        for i=1:size(conedges,2)
            t1 = find(t(1,:) == conedges(1,i));
            if isempty(t1)
                t1 = 0;
            end
            t2 = find(t(2,:) == conedges(1,i));
            if isempty(t2)
                t2 = 0;
            end
            t3 = find(t(3,:) == conedges(1,i));
            if isempty(t3)
                t3 = 0;
            end
            %aux = intersect(t1,t2);
            %n1 = intersect(aux,t3)
            n1 = [t1,t2,t3];

            t1 = find(t(1,:) == conedges(2,i));
            if isempty(t1)
                t1 = 0;
            end
            t2 = find(t(2,:) == conedges(2,i));
            if isempty(t2)
                t2 = 0;
            end
            t3 = find(t(3,:) == conedges(2,i));
            if isempty(t3)
                t3 = 0;
            end
            %aux = intersect(t1,t2);
            %n2 = intersect(aux,t3)

            n2 = [t1,t2,t3];
            aux = intersect(n1,n2);
            ind = aux>0 ;
            triangule(i) = aux(ind);
        end
    end
    
    
    
    
    

end