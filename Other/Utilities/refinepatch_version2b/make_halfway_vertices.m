function [Vout,HT_index, HT_values]=make_halfway_vertices(EV_table,ET_table,ETV_index,V,Ne)
% ETV_index_sparse=sparse(size(V,1),size(V,1));
%
% for i=1:size(ETV_index,1)
% 	if(ETV_index(i,1)>0)
% 		ETV_index_sparse(ETV_index(i,1),ETV_index(i,2))=i;
% 	else
% 		ET_table=ET_table(1:(i-1),:);
% 		EV_table=EV_table(1:(i-1));
% 		break
% 	end
% end


% Table to cell array
ETV_index_vall=cell(length(V),1);
for i=1:size(ETV_index,1)
    if(ETV_index(i,1)>0)
        ETV_index_vall{ETV_index(i,1)}=[ETV_index_vall{ETV_index(i,1)} i];
    else
        ET_table=ET_table(1:(i-1),:);
        EV_table=EV_table(1:(i-1));
        break;
    end
end

HT_index=cell(length(V),1);
HT_values=cell(length(V),1);

% Make output V
Vout=zeros(size(V,1)*4,3);
Vout(1:size(V,1),:)=V;
Vindex=size(V,1);

for i=1:length(V)
    Pneig=Ne{i};
    for j=1:length(Pneig);
        % Get the tangent and velocity of the edge P -> Pneig
        index=Ne{i}; vals=ETV_index_vall{i}; select1=vals(index==Pneig(j));
        Va=EV_table( select1); Ea=ET_table( select1,:);
        % Get the tangent and velocity of the edge Pneig -> P
        index=Ne{Pneig(j)}; vals=ETV_index_vall{Pneig(j)}; select2=vals(index==i);
        Vb=EV_table(select2); Eb=ET_table(select2,:);
            
        % The four points describing the spline
        P0=V(i,:);
        P3=V(Pneig(j),:);
        P1=P0+Ea*Va/3;
        P2=P3+Eb*Vb/3;
               
        % Spline used to calculated the xyz coordinate of the middle of each edge;
        c = 3*(P1 - P0);
        b = 3*(P2 - P1) - c;
        a = P3 - P0 - c - b;
       
        halfwayp = a*0.125 + b*0.250 + c*0.500 + P0;
        
        % Save the edge middle point
        if(sum(HT_index{i}==Pneig(j))==0)
            Vindex=Vindex+1;
            Vout(Vindex,:)=halfwayp;
            HT_index {i}= [HT_index{i} Pneig(j)];
            HT_values{i}=[HT_values{i} Vindex];
            HT_index {Pneig(j)}=[HT_index{ Pneig(j)} i];
            HT_values{Pneig(j)}=[HT_values{Pneig(j)} Vindex];
        end
    end
end
Vout=Vout(1:Vindex,:);
