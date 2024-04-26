function [V2 F2 ID2] = InsidePoints3D(V,F)
% Remove points outside the contour
pointOutside=true(size(V,1),1);
for i=1:size(F,1)
    pointOutside(F(i,1))=false; 
    pointOutside(F(i,2))=false;
    pointOutside(F(i,3))=false;
end

for i=1:size(V,1)
    if(pointOutside(i))
        if(CheckPointOutInside3D(V(i,:),V,F)), pointOutside(i)=false; end
    end
end
V2=[V (1:size(V,1))'];
V2(pointOutside,:)=[];
ID2=V2(:,4);  V2=V2(:,1:3); 

F2=F; for j=1:length(V2), F2(F==ID2(j))=j; end
