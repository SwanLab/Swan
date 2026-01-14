function  [cX,bknn]= getcX(stX,k,hknn)                     
cluster=[]; n=length(stX);B=1; i=1;           
 knn=1;bknn=[];cX=[];ccY=[];
while(knn<hknn) %hknn是knn的上界。
     options = struct( 'NN', knn, 'GraphDistanceFunction', 'euclidean');                  
    W = adjacency(options,stX);
    a=getjuzhen(W);
    [cluster]=getcu(a); 
    cluster(cellfun('isempty',cluster)) = [];
    K=length(cluster);
    if (K>k)
       knn=knn+2;
       continue;
    else if (K==k)
        bknn=knn;
        else 
            break;
        end
    end   
     for i=1:k
       Xi=stX(cluster{1,i},:);
       ccY(1,1: length(cluster{1,i}))=i;
       ccX=([ccY' Xi]); 
       ccY=[];
       cX=([cX ;ccX]); 
     end
     break;
   bknn;
end
if isempty(bknn)
    knn=1;  options = struct( 'NN', knn, 'GraphDistanceFunction', 'euclidean');                  
    W = adjacency(options,stX);
    a=getjuzhen(W);
    [cluster]=getcu(a);  cluster(cellfun('isempty',cluster)) = []; K=length(cluster); 
           Mi=zeros(K,K);dX=[];
        for j=1:K
              Xj=stX(cluster{1,j},:);
              mj=mean(Xj);
              dX=[dX mj];
            for i=(j+1):K
                    Xi=stX(cluster{1,i},:);
                    mi=mean(Xi);
                   Mi(j,i)=sum((mj-mi).^2,2);
            end
        end
                [r,c]=find(Mi==min(nonzeros(triu(Mi,1))),1);
                cluster{1,r}=([cluster{1,r} cluster{1,c}]);
                cluster{1,c}=[];
                 b=cellfun('isempty',cluster);  
                K1=length( find( b(1,:)==0 ) );
             while(K1~=k)   
               Mi(c,:)=NaN;  Mi(:,c)=NaN;dX(:,c)=NaN;
               Xj1=stX(cluster{1,r},:);
               mr=mean(Xj1);
               for p=1:K
                   Mi(r,p)=sum((mr-dX(:,p)).^2,2); 
               end
                Mi(r,c)=NaN;
                [r,c]=find(Mi==min(min(Mi)),1);
                cluster{1,r}=([cluster{1,r} cluster{1,c}]);
                cluster{1,c}=[];
                b=cellfun('isempty',cluster);  
                K1=length( find( b(1,:)==0 ) );
             end
           cluster(cellfun('isempty',cluster)) = [];     
    for i=1:k
          Xi=stX(cluster{1,i},:);
          ccY(1,1: length(cluster{1,i}))=i;
          ccX=([ccY' Xi]); 
          ccY=[];
          cX=([cX ;ccX]); 
    end
end 
end
function a=getjuzhen(W)
[m,~]=size(W);
a=zeros(m,m);
[r,c]=find((W+W')>0);
for i=1:length(r)
a(r(i,1),c(i,1))=1;
end
end
