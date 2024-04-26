function [V,F,nF,T,nT]=process(V,F,nF,T,nT,Frows,F_Local,F_Local_new, Vertex_id)

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

    nT=nT+1;
    Tc=delaunayn(V([F_Local_new(i,:) Vertex_id],:));
    Tc2=Tc; 
    Tc2(Tc==1)=F_Local_new(i,1); 
    Tc2(Tc==2)=F_Local_new(i,2); 
    Tc2(Tc==3)=F_Local_new(i,3); 
    Tc2(Tc==4)=Vertex_id;
    T(nT,:)=Tc2; 
end
