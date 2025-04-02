function [ET_table,EV_table,ETV_index]=edge_tangents_double(V,Ne)
% Edge tangents table
ET_table=zeros(size(V,1)*4,3);
% Edge velocity table
EV_table=zeros(size(V,1)*4,1);
% Edge tangents/velocity index for tables
ETV_index=zeros(size(V,1)*4,2);
ETV_num=0;

% Calculate the tangents and velocity for each edge
Pn=zeros(length(Ne),3); Pnop=zeros(length(Ne),3);
for i=1:size(V,1) 
    P=V(i,:);
    Pneig=Ne{i};

    % Find the opposite vertex of each neigbourh vertex.
    % incase of odd number of neigbourhs interpolate the opposite neigbourh
    if(mod(length(Pneig),2)==0)
       for k=1:length(Pneig)
            neg=k+length(Pneig)/2; if(neg>length(Pneig)), neg=neg-length(Pneig); end
            Pn(k,:) = V(Pneig(k),:); Pnop(k,:) = V(Pneig(neg),:);
       end
    else
       for k=1:length(Pneig)
           neg=k+length(Pneig)/2; neg1=floor(neg); neg2=ceil(neg);
           if(neg1>length(Pneig)), neg1=neg1-length(Pneig); end
           if(neg2>length(Pneig)), neg2=neg2-length(Pneig); end
           Pn(k,:) = V(Pneig(k),:); Pnop(k,:) = (V(Pneig(neg1),:)+V(Pneig(neg2),:))/2;
       end
    end

    for j=1:length(Pneig);
        % Calculate length edges of face
        Ec=sqrt(sum((Pn(j,:)-P).^2))+1e-14;  
        Eb=sqrt(sum((Pnop(j,:)-P).^2))+1e-14; 
        Ea=sqrt(sum((Pn(j,:)-Pnop(j,:)).^2))+1e-14; 

        % Calculate face surface area
        s = ((Ea+Eb+Ec)/2); 
        h = (2/Ea)*sqrt(s*(s-Ea)*(s-Eb)*(s-Ec))+1e-14;
        x = (Ea^2-Eb^2+Ec^2)/(2*Ea);
 
        % 2D triangle coordinates
        % corx(1)=0;    cory(1)=0;
        % corx(2)=x;    cory(2)=h;
        % corx(3)=Ea;   cory(3)=0;
        % corx(4)=0;    cory(4)=0;
       
        % Calculate tangent of 2D triangle
        Np=[-h x]; Np=Np/(sqrt(sum(Np.^2))+1e-14); 
        Ns=[h Ea-x]; Ns=Ns/(sqrt(sum(Ns.^2))+1e-14); 
        Nb=Np+Ns; 
        Tb=[Nb(2) -Nb(1)];
        
        % Back to 3D coordinates
        Pm=(Pn(j,:)*x+Pnop(j,:)*(Ea-x))/Ea;
        X3=(Pn(j,:)-Pnop(j,:))/Ea;
        Y3=(P-Pm)/h;
   
        % 2D tangent to 3D tangent
        Tb3D=(X3*Tb(1)+Y3*Tb(2));  Tb3D=Tb3D/(sqrt(sum(Tb3D.^2))+1e-14); 
        
        % Edge Velocity
        Vv=0.5*(Ec+0.5*Ea);
        
        ETV_num=ETV_num+1;
        ETV_index(ETV_num,:)=[i Pneig(j)];
        ET_table(ETV_num,:)= Tb3D;
        EV_table(ETV_num)=Vv;
    end
end

