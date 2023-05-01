elementNumberX = 200;
elementNumberY = 100;
nodesNumerationVec = 1:(elementNumberX+1)*(elementNumberY+1);
nodesNumerationMat = reshape(nodesNumerationVec,(elementNumberY+1),(elementNumberX+1));
[y,x] = size(nodesNumerationMat);
connectivityMatrix  = zeros(elementNumberX*elementNumberY,6);
cont = 1;
for xE = 1:x-1
    for yE = 1:y-1
         connectivityMatrix(cont,1) = cont;
         connectivityMatrix(cont,2) = nodesNumerationMat(yE+1,xE);
         connectivityMatrix(cont,3) = nodesNumerationMat(yE+1,xE+1);
         connectivityMatrix(cont,4) = nodesNumerationMat(yE,xE+1);
         connectivityMatrix(cont,5) = nodesNumerationMat(yE,xE);
         cont = cont+1;
    end 
end 
nodesCoords = zeros(x*y,4);
cont = 1;
xCoord = 0;
yCoord = (y)/100;
for xE = 1:x
    for yE = 1:y
         nodesCoords(cont,1) = cont;
         nodesCoords(cont,2) = xCoord;
         nodesCoords(cont,3) = yCoord;
         cont = cont+1;
         yCoord = yCoord-0.01;
    end 
    yCoord =(elementNumberY+1)/100;
    xCoord = xCoord+0.01;
end 
nodesCoords(nodesCoords <0) = 0;

prescribedNode = zeros(y+2,3);
for yE = 1:y
    prescribedNode(yE,1) = yE;
    prescribedNode(yE,2) = 1;
    prescribedNode(yE,3) = 0;
end
prescribedNode(y+1,:)=[(x)*(y) 1 0];
prescribedNode(y+2,:)=[(x)*(y) 2 0];
imposedForce = [1 1 1e-3];

pointload_complete = imposedForce;
lnodes = prescribedNode;
gidcoord = nodesCoords;
gidlnods = connectivityMatrix; 

save('imposedForce.mat','pointload_complete');
save('prescribedNode.mat','lnodes');
save('coordsMatrix.mat','gidcoord');
save('connectivityMatrix.mat','gidlnods');
% fileID = fopen('connectivityMatrix.txt','w');
% fprintf(fileID,'%d %d %d %d %d %d\n',connectivityMatrix');
% fclose(fileID);
% fileID = fopen('coordsMatrix.txt','w');
% fprintf(fileID,'%d %d %d %d\n',connectivityMatrix');
% fclose(fileID);

