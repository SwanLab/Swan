function volume=CheckVolumeTetraMesh(V,T)
volume=0;
for i=1:size(T,1)
    a=V(T(i,1),:);
    b=V(T(i,2),:);
    c=V(T(i,3),:);
    d=V(T(i,4),:);

    a=a-d;
    b=b-d;
    c=c-d;

    k=cross(b,c);
    v = abs(a(1)*k(1)+a(2)*k(2)+a(3)*k(3))/6;

    volume=volume+v;
end


function c=cross(a,b)
a=a(:); b=b(:); 
c = [a(2,:).*b(3,:)-a(3,:).*b(2,:)
     a(3,:).*b(1,:)-a(1,:).*b(3,:)
     a(1,:).*b(2,:)-a(2,:).*b(1,:)];
 

