function PlotMesh(T,X,str,nonum)
% PlotMesh(T,X,str,nonum)
% X:    nodal coordinates
% T:    connectivities 
% str:  linestyle, color and marker used in the plot (optional)
% nonum = 1 to show nodes' number(optional)
%


% Line style and color
if nargin == 2
  str1 = 'yo';
  str2 = 'y-';
else
  if str(1) == ':' | str(1) == '-'  
    str1 = 'yo'; 
    str2 = ['y' str];
  else
    str1 = [str(1) 'o'];
    str2 = str;
  end
end
 
nnodes = size(T,2);


if  size(X,2) == 2
    order = [1:nnodes,1];
    if nnodes == 8
         order = [1 5 2 6 3 7 4 8 1];
    end;
    % Nodes
    plot(X(:,1),X(:,2),str1)
    hold on
    % Elements
    for j = 1:size(T,1)
        plot(X(T(j,order),1),X(T(j,order),2),str2)
    end
end

% nodes number
if nargin==4 
    if nonum==1
        for I=1:size(X,1)
            text(X(I,1),X(I,2),int2str(I))
        end
    end
end

axis('equal')    
axis('off') 

hold off 
