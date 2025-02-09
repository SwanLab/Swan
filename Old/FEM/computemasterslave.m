
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
    for j=1:size(minnodes)
        cond1=gidcoord(maxnodes,idime_next)==gidcoord(minnodes(j),idime_next);
        cond2=gidcoord(maxnodes,idime_next2)==gidcoord(minnodes(j),idime_next2);
         k=cond1.*cond2 ==1 ;
         Master_slave=[Master_slave;[minnodes(j),maxnodes(k)]];
    end    
end
i=0;
j=0;
k=0;
for inode=1:size(gidcoord,1)
    if gidcoord(inode,3)==0 && gidcoord(inode,4)==0
        i=i+1;
        xaxis(i)=inode;
    end
    if gidcoord(inode,2)==0 && gidcoord(inode,4)==0
        j=j+1;
        yaxis(j)=inode;
    end
    if gidcoord(inode,2)==0 && gidcoord(inode,3)==0
        k=k+1;
        zaxis(k)=inode;
    end
end
for ix=1:length(xaxis)
   check1=gidcoord(:,2:4)==[gidcoord(xaxis(ix),2),0,1];
   check2=gidcoord(:,2:4)==[gidcoord(xaxis(ix),2),1,1];
   check3=gidcoord(:,2:4)==[gidcoord(xaxis(ix),2),1,0];
   check1=sum(check1,2)==3;
   check2=sum(check2,2)==3;
   check3=sum(check3,2)==3;
   matrix=[repmat(xaxis(ix),3,1),[gidcoord(check1,1);gidcoord(check2,1);gidcoord(check3,1)]];
   null_elements=any(ismember(Master_slave,matrix)');
   Master_slave(null_elements,:)=[];
   Master_slave=[Master_slave;matrix];  
end
for ix=1:length(yaxis)
   check1=gidcoord(:,2:4)==[0,gidcoord(yaxis(ix),3),1];
   check2=gidcoord(:,2:4)==[1,gidcoord(yaxis(ix),3),1];
   check3=gidcoord(:,2:4)==[1,gidcoord(yaxis(ix),3),0];
   check1=sum(check1,2)==3;
   check2=sum(check2,2)==3;
   check3=sum(check3,2)==3;
   matrix=[repmat(yaxis(ix),3,1),[gidcoord(check1,1);gidcoord(check2,1);gidcoord(check3,1)]];
   null_elements=any(ismember(Master_slave,matrix)');
   Master_slave(null_elements,:)=[];
   Master_slave=[Master_slave;matrix];  
end
for ix=1:length(zaxis)
   check1=gidcoord(:,2:4)==[0,1,gidcoord(zaxis(ix),4)];
   check2=gidcoord(:,2:4)==[1,1,gidcoord(zaxis(ix),4)];
   check3=gidcoord(:,2:4)==[1,0,gidcoord(zaxis(ix),4)];
   check1=sum(check1,2)==3;
   check2=sum(check2,2)==3;
   check3=sum(check3,2)==3;
   matrix=[repmat(zaxis(ix),3,1),[gidcoord(check1,1);gidcoord(check2,1);gidcoord(check3,1)]];
   null_elements=any(ismember(Master_slave,matrix)');
   Master_slave(null_elements,:)=[];
   Master_slave=[Master_slave;matrix];  
end
null_elements=any(ismember(Master_slave,lnodes(:,1))');
Master_slave(null_elements,:)=[];