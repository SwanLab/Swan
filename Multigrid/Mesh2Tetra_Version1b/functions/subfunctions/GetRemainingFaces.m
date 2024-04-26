function Faces=GetRemainingFaces(T,F,V)

Faces1=T(:,[3 2 1]);
Faces2=T(:,[2 4 1]);
Faces3=T(:,[4 3 1]);
Faces4=T(:,[4 2 3]);

Faces=[Faces1;Faces2;Faces3;Faces4;F];
FacesS=sort(Faces,2);

Faces=[FacesS(:,1)*(size(V,1)^2)+FacesS(:,2)*size(V,1)+FacesS(:,3) Faces];
Faces=sortrows(Faces,1);

FacesOld=0;
Facesl=false(size(Faces,1),1);
for i=1:size(Faces,1)
    if(Faces(i,1)==FacesOld)
        Facesl(i-1)=true; 
        Facesl(i)=true; 
    end
    FacesOld=Faces(i,1);
end
Faces(Facesl,:)=[];
Faces=Faces(:,2:4);
