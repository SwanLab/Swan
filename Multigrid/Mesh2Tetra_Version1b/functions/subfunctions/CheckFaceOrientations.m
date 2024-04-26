function inter=CheckFaceOrientations(F)
Edges1=F(:,[1 2]); Edges2=F(:,[2 3]); Edges3=F(:,[3 1]);
Edges=[Edges1;Edges2;Edges3];
Edges2=unique(Edges,'rows');

inter=false;
for i=1:size(Edges2,1)
    Nn=sum((Edges2(i,1)==Edges(:,1))&(Edges2(i,2)==Edges(:,2)));
    Ni=sum((Edges2(i,1)==Edges(:,2))&(Edges2(i,2)==Edges(:,1)));
    if(Nn~=Ni), inter=true; break; end
end

