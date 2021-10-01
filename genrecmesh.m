function [conectividades,nodes] = genrecmesh(xi,yi)
%------------------------------ Lattice, truss/frame mesh -------------------------------%
% Made by : Diego Alfredo Quexada Rodríguez
% Colaboration with Jair Tovar.
% Universidad Nacional de Colombia
% This script generates a rectangular mesh of size x by y, where nodes and
% conectivities are exported to a text file for further analysis with FEM.
% Elements are seen as lines wich are arranged in squares with
% diagonals between its vertices. The files generated are to be read in
% abaqus. Hope it helps.
%-----------------------------------------------------------------------C--%
nelemx=size(xi,2)-1;
nelemy=size(yi,2)-1;
elemh=zeros(nelemx*(nelemy+1),3);
elemv=zeros((nelemx+1)*nelemy,3);
elemd1= zeros(nelemx*nelemy,3); 
elemd2= zeros(nelemx*nelemy,3); 
cont=1;
%CICLO HORIZONTALES
for i=1:nelemy+1
    for j=1:nelemx
        elemh(cont,1)= cont;
        elemh(cont,2)= j+(nelemx+1)*(i-1);
        elemh(cont,3)= j+(nelemx+1)*(i-1)+1;
        cont=cont+1;  
    end
end
%CICLO VERTICALES
cont2=1;
for i=1:nelemy
    for j=1:nelemx+1
        
        elemv(cont2,1)= cont;
        elemv(cont2,2)= j+(nelemx+1)*(i-1);
        elemv(cont2,3)= j+(nelemx+1)+(nelemx+1)*(i-1);
        cont2=cont2+1;
      cont=cont+1;
    end
end
% %diagonal cycle 1 ///
cont2=1;
for i=1:nelemy
    for j=1:nelemx
        
        elemd1(cont2,1)= cont;
        elemd1(cont2,2)= j+(nelemx+1)*(i-1);
        elemd1(cont2,3)= j+(nelemx+1)+(nelemx+1)*(i-1)+1;
        cont2=cont2+1;
      cont=cont+1;
    end
end
% diagonal cycle 2  \\\\
cont2=1;
for i=1:nelemy
    for j=1:nelemx
        
        elemd2(cont2,1)= cont;
        elemd2(cont2,2)= j+(nelemx+1)*(i-1)+1;
        elemd2(cont2,3)= j+(nelemx+1)+(nelemx+1)*(i-1);
        cont2=cont2+1;
      cont=cont+1;
    end
end
conectividades=[elemh;elemv;elemd1;elemd2];
numnodes=length(xi)*length(yi);
nodes=zeros(numnodes,1);
zeroz=zeros(numnodes,1);
for i=1:numnodes
    nodes(i,1)=i;
end
conectividades = conectividades(all(conectividades,2),:);%deleting 0 rows
for i=1:length(conectividades(:,1))
    conectividades(i,1)=i;
end
   %Creating node matrix
[A,B] = meshgrid(xi,yi);
c=cat(2,A',B');
d=reshape(c,[],2);
nodes=[nodes d];
