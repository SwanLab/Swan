function [V,F,nF,T,nT]=collapse_edge(V,F,nF,T,nT,Volume_Original,Options)
% Make list with all vertices which have not collapsed yet
Vind=make_left_vertice_list(V,F,nF,Options) ;
collapse_possible=false;
vs=0;
while(~collapse_possible)
    % Select the first vertice
    vs=vs+1;
    if(vs>length(Vind)), return, end
    Vertex_id=Vind(vs);
    
    % Find the faces involved
    Irows=F(1:nF,:)==Vertex_id;
    [Frows,~] = find(Irows);
    
    % Local Patch
    F_Local=F(Frows,:);
    
    % Find Local neighbour vertices
    Vind_Local=unique(F_Local); 
    [~,ind]=sort(rand(size(Vind_Local)));
   % Vind_Local=Vind_Local(ind);

    Vind_Local(Vind_Local==Vertex_id)=[];
    
    for il=1:size(Vind_Local)
        % Select a vertice from the Local patch list
        Vertex_id_Local=Vind_Local(il);
        
        % Do the edge collapse
        F_Local_new=F_Local; F_Local_new(F_Local_new==Vertex_id)=Vertex_id_Local;
        
        % Remove the faces collapsed to an edge
        col_edge=(F_Local_new(:,1)==Vertex_id_Local)+(F_Local_new(:,2)==Vertex_id_Local)+(F_Local_new(:,3)==Vertex_id_Local);
        F_Local_new(col_edge>1,:)=[];
        
        % Check if the new edge is inside
        check(1)=(size(F_Local_new,1)==0);
        if(check(1)),continue; end
        [inter,F2,nF2,Tnew]=VolumeCheckNew(V,F,nF,T,nT,Frows,F_Local,F_Local_new, Vertex_id,Vertex_id_Local,Volume_Original);
        check(2)=inter;
        if(check(2)),continue; end
        check(3)=CheckFaceOrientations(F2(1:nF2,:));
        if(check(3)), continue; end
        check(4)=~CheckMoveInside3D(V,F_Local_new,Vertex_id);
        if(check(4)), continue; end
        check(5)=CheckMeshInterSections(V,F2(1:nF2,:)); 
        if(check(5)), continue; end
        %check(6)=CheckPointsInsideTetrahedron(V,F_Local_new, Vertex_id);
        %if(check(6)), continue; end
        if(~isempty(Tnew))
            F1=Tnew(:,[1 2 3]); F2=Tnew(:,[1 4 2]); F3=Tnew(:,[1 3 4]); F4=Tnew(:,[3 2 4]);  F5=F(1:nF,:);
            Ft=[F1;F2;F3;F4;F5];
            check(7)=CheckMeshInterSections(V,Ft,size(Tnew,1)*4); 
            if(check(7)), continue; end
        end
        
        check=~any([check(1) check(2) check(3) check(4) check(5) check(6)]);
        if(check), collapse_possible=true; break; end
    end
end
[V,F,nF,T,nT]=process(V,F,nF,T,nT,Frows,F_Local,F_Local_new, Vertex_id);