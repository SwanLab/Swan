function Visible=CheckPointOutInside3D(P3,V,F)
ka=0; kb=0; kc=0; lista=[]; listb=[]; listc=[];
P4a=P3+[ 123456781 9831542342 5831542342];
P4b=P3+[-124342343 3234454234 3831542342];
P4c=P3+[ 648832349 5415435923 -1831542342];

for k=1:size(F,1)
    CF=F(k,:); 
    O1=V(CF(1),:); O2=V(CF(2),:);  O3=V(CF(3),:);
        
    [inter,inter_x,inter_y,inter_z]=LineTriangleIntersection(O1,O2,O3,P3,P4a);
    if(inter), ka=ka+1; lista(ka,:)=[inter_x inter_y inter_z]; end

    [inter,inter_x,inter_y,inter_z]=LineTriangleIntersection(O1,O2,O3,P3,P4b);
    if(inter), kb=kb+1; listb(kb,:)=[inter_x inter_y inter_z]; end
    
    [inter,inter_x,inter_y,inter_z]=LineTriangleIntersection(O1,O2,O3,P3,P4c);
    if(inter), kc=kc+1; listc(kc,:)=[inter_x inter_y inter_z]; end
end

lista=unique(round(lista*1e8)/1e8,'rows');
listb=unique(round(listb*1e8)/1e8,'rows');
listc=unique(round(listc*1e8)/1e8,'rows');
check=[mod(size(lista,1),2)==0 mod(size(listb,1),2)==0 mod(size(listc,1),2)==0];
if(sum(check)>1),
    Visible=false; 
else
    Visible=true; 
end
