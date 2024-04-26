function Visible=CheckVisiblePoint3D(V,F,i,j)
Visible=true;
P1=V(i,:); 
P2=V(j,:);


if(any(any(F==i,2)&any(F==j,2))), return, end


for k=1:size(F,1)
    CF=F(k,:);
    if(any(CF==i)||any(CF==j)), continue, end
    O1=V(CF(1),:); O2=V(CF(2),:);  O3=V(CF(3),:);
    
    inter=LineTriangleIntersection(O1,O2,O3,P1,P2,true);

     
%     FV.vertices=V;
%     FV.faces=CF;
%     subplot(3,4,k),hold on; patch(FV,'facecolor',[1 0 0]);
%     plot3([V(i,1) V(j,1)],[V(i,2) V(j,2)],[V(i,3) V(j,3)]);

    if(inter), Visible=false; return; end
end


check=any(any(F==i,2)&any(F==j,2));
if(~check)
    P3=(V(i,:)+V(j,:))/2;
    Visible=CheckPointOutInside3D(P3,V,F);
end
