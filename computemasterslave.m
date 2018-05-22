
Master_slave=[];
for idime=2:4
    maxnodes=find(gidcoord(:,idime)==max(gidcoord(:,idime)));
    minnodes=find(gidcoord(:,idime)==min(gidcoord(:,idime)));    
    if idime==2
        idime_next=3;
        idime_next2=4;
    end    
    if idime==3
        idime_next=4;
        idime_next2=2;
    end
    if idime==4
        idime_next=2;
        idime_next2=3;
    end  


    for j=1:size(maxnodes)
        cond1=gidcoord(minnodes,idime_next)==gidcoord(maxnodes(j),idime_next);
        cond2=gidcoord(minnodes,idime_next2)==gidcoord(maxnodes(j),idime_next2);
         k=cond1.*cond2 ==1 ;
         Master_slave=[Master_slave;[maxnodes(j),minnodes(k)]];
    end
    
end

null_elements=any(ismember(Master_slave,lnodes(:,1))');
Master_slave(null_elements,:)=[];