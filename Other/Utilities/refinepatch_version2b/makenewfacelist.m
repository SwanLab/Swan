function Fnew=makenewfacelist(F,HT_index, HT_values)
% Combine the edge middle points and old vertex points to faces.
% (4 Split method)
Fnew=zeros(length(F)*4,3);
for i=0:length(F)-1,
    vert1=F(i+1,1); 
    vert2=F(i+1,2); 
    vert3=F(i+1,3);
    
    index=HT_index{vert1}; vals=HT_values{vert1};
    verta= vals(index==vert2);
    index=HT_index{vert2}; vals=HT_values{vert2};
    vertb= vals(index==vert3);
    index=HT_index{vert3}; vals=HT_values{vert3};
    vertc= vals(index==vert1);

    Fnew(i*4+1,:)=[vert1 verta vertc];
    Fnew(i*4+2,:)=[verta vert2 vertb];
    Fnew(i*4+3,:)=[vertc vertb vert3];
    Fnew(i*4+4,:)=[verta vertb vertc];
end
