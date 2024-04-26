function volume=CheckVolumeFaceMesh(V,F)
volume=0;
for i=1:size(F,1)
    a=V(F(i,1),:); b=V(F(i,2),:); c=V(F(i,3),:);

    k=cross(b,c);
    v = (a(1)*k(1)+a(2)*k(2)+a(3)*k(3))/6;
    
    volume=volume+v;
end
volume=-(volume);

function c=cross(a,b)
a=a(:); b=b(:); 
c = [a(2,:).*b(3,:)-a(3,:).*b(2,:)
     a(3,:).*b(1,:)-a(1,:).*b(3,:)
     a(1,:).*b(2,:)-a(2,:).*b(1,:)];
 