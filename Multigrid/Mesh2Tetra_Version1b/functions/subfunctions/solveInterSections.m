function [F,V]=solveInterSections(F,V)
FV.faces=F; FV.vertices=V;
nInter=inf; nInterOld=inf;
while(nInter>0)
    [FV,nInter]=solveInterSectionsSub(FV);
    if((nInter>0)&&(nInter==nInterOld))
        error('solveInterSections:process','Mesh intersections not solved');
    end
    nInterOld=nInter;
end
F=FV.faces; V=FV.vertices;


function [FV,nInter]=solveInterSectionsSub(FV)
[intersect FV2]=CheckMeshInterSections(FV.vertices,FV.faces,[],true);
if(intersect), nInter=size(FV2.faces,1); else nInter=0; return; end
    
Vind=unique(FV2.faces);
for i=1:length(Vind);
    Vertex_id=Vind(i);

    % Find the faces involved
    Irows=FV.faces==Vertex_id;
    [Frows,~] = find(Irows);

    % Local Patch
    F_Local=FV.faces(Frows,:);

    % Find Local neighbour vertices
    Vind_Local=unique(F_Local); 
    Vind_Local(Vind_Local==Vertex_id)=[];

    for il=1:size(Vind_Local)
        [nInter_new FV_new]=trySollution(il,Vind_Local,F_Local,Frows,FV,Vertex_id);
        if(nInter_new<nInter)
            FV=FV_new;
            return
        end
    end
end

function [nInter FV]=trySollution(il,Vind_Local,F_Local,Frows,FV,Vertex_id)
    F=FV.faces;
    nF=size(F,1);
    
    % Select a vertice from the Local patch list
    Vertex_id_Local=Vind_Local(il);

    % Do the edge collapse
    F_Local_new=F_Local; F_Local_new(F_Local_new==Vertex_id)=Vertex_id_Local;

    % Remove the faces collapsed to an edge
    col_edge=(F_Local_new(:,1)==Vertex_id_Local)+(F_Local_new(:,2)==Vertex_id_Local)+(F_Local_new(:,3)==Vertex_id_Local);
    F_Local_new(col_edge>1,:)=[];

    % Remove faces by replacing them
    Frows=sort(Frows(:),1,'descend');
    for i=1:size(F_Local,1),
        F(Frows(i),:)=F(nF,:); nF=nF-1;
    end

    for i=1:size(F_Local_new,1)
        F2=F(1:nF,:);
        checkduplicate=(any(F2==F_Local_new(i,1),2)&any(F2==F_Local_new(i,2),2)&any(F2==F_Local_new(i,3),2));
        if(any(checkduplicate))
            j=find(checkduplicate);
            if(length(j)>1)
                error('boundary_collapse:dublicates','This is not possible');
            end
            F(j,:)=F(nF,:); nF=nF-1;
        else
            nF=nF+1; F(nF,:)=F_Local_new(i,:);
        end
    end

    FV.faces=F(1:nF,:);
    [intersect FV3]=CheckMeshInterSections(FV.vertices,FV.faces,[],true);
    if(intersect),nInter=size(FV3.faces,1); else nInter=0; end