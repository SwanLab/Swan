% LkPPC
% % A Matlab code for Local k-Proximal Plane Clustering
% Reference
% % Yuan-Hai Shao, Yan-Ru Guo, 
% Main Function
% % Need stdata function; adjacency function; GepOneSide function;
% % getcu function ; getcX function ; GetchushiW function;
% %%% function pY =LkPPC(stX,cX,k,c,g,W)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % LocalPPC: Local Proximal Plane Clustering
% % % 
% % %pY =LkPPC(stX,cX,k,c,g,W)
% % % Input:  stX -  data points matrix .Each row vector of fea is a data point and is normalized with the mean 0 and standard deviation 1. 
% % %         cX-  the data points who is used to construct the initial plane
% % %          k-   number of cluster;
% % %          hknn- the upper bound of the KNN;
% % %           W- the construct the initial plane;
% % %         FunPara - Struct value in Matlab. The fields in options that can be set:
% % %              c: [0,inf] appropriate parameter to tune the weight. 
% % %              g:[0,1]  is used to control the localization of the clustering plane.
% % %              Output:
% % % Output:  pY - Predict the class of X.
% Examples
   X = rand(50,10);Y=[ones(20,1); ones(20,1)+1; ones(10,1)+2];
    c=0.01; g=0.1; k=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initailization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[s,t]=size(X);hknn=20;B=1;knn=1; bknn=[];cX=[];ccY=[];i=1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each data point is normalized with the mean 0 and standard deviation 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 stX=stdata(X);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use KNN to find cX-  the data points who is used to construct the initial plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [cX,bknn]= getcX(stX,k,hknn);    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construct the initial plane W
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  W=GetchushiW(cX,k,c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update plane W:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pYnew=zeros(s,1);[s,t]=size(stX);[mm,nn]=size(cX);pY=crossvalind('kfold',s,k);
z=0;q=zeros(k,t); V=zeros(s,k);
while(~isempty(find(pY~=pYnew, 1)) && z~=1000)   
    pYnew=pY;
    z=z+1;%控制运行的次数。
    % update W
     for i=1:k
             tA=stX((pY==i),:);
             tB=stX((pY~=i),:); 
             mi=size(tA,1);
         if ~isempty(find(pY==i, 1))
            W(i,:)=GepOneSide(tA,tB,c); 
             q(i,:)=sum(tA)/mi;  
         end         
     end   
    for l=1:s
         for ff=1:k    
             V(l,ff)=(norm((stX(l,:)-q(ff,:)),2)^2);
         end
     end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predict and output pY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
     pY=abs(stX*W(:,1:t)'+ones(s,1)*W(:,t+1)')+g*V(:,:);
               [tmp,pY]=min(pY');
                pY=pY';
end

    

         