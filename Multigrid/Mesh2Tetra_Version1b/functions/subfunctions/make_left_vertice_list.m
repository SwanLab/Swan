function Vind=make_left_vertice_list(V,F,nF,Options)
% Make list with all vertices which have not collapsed yet
Vind=unique(F(1:nF,:));

%sortvalue=zeros(size(Vind));
% for i=1:length(Vind)
%     % Get Neighbours
%     ne=unique(F(any(F(1:nF,:)==Vind(i),2),:)); 
%     Vne= V(ne,:);
%     
%     dx=repmat(Vne(:,1),1,size(Vne,1))-repmat(Vne(:,1)',size(Vne,1),1);
%     dy=repmat(Vne(:,2),1,size(Vne,1))-repmat(Vne(:,2)',size(Vne,1),1);
%     dz=repmat(Vne(:,3),1,size(Vne,1))-repmat(Vne(:,3)',size(Vne,1),1);
%     d=sqrt(dx.^2+dy.^2+dz.^2);
%     dc=d(ne==Vind(i),:);
%     dn=mean(d(:));
%     dc=mean(dc(:));
%     sortvalue(i)=dn/(dc+eps);
% end
if(Options.mode==1)
    sortvalue=rand(size(Vind));
    [~,ind]=sort(sortvalue);
    Vind=Vind(ind);
end


