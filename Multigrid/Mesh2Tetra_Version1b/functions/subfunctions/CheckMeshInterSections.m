function [intersect FV]=CheckMeshInterSections(V,F,nMax,findall)
FV=[];
nF=size(F,1);
if((nargin<3)||isempty(nMax))
    nMax=nF-1;
else
    nMax=min(nMax,nF-1);
end
if(nargin<4)
    findall=false; 
end

FV.vertices=V;
intersect=false;
k=1;
for j=1:nMax
    O1=V(F(j,1),:); O2=V(F(j,2),:); O3=V(F(j,3),:); 
    for i=j+1:nF
        P1=V(F(i,1),:); P2=V(F(i,2),:); P3=V(F(i,3),:); 
        if(all(F(i,:)==F(j,:))), continue; end
        
        inter=TriangleTriangleIntersection(P1,P2,P3,O1,O2,O3,true);
         if(inter), 
             FV.faces(k,:)=F(j,:); k=k+1;
             FV.faces(k,:)=F(i,:); k=k+1;
             intersect=true; 
             if(~findall)
                break;
             end
         end
    end
    if(intersect), break; end
end

