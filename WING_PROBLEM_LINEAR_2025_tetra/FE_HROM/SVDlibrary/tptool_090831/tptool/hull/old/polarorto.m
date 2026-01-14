%az U1 ponthalmazt burkolja polarv alaku tetraéderrel
%fi2 a tetraéder csúcsai
% U1=U2*fi2

function [U2,v]=polarorto(U1,polarv)

[n r]=size(U1);

a=polarv;

for i=1:r-2    
    v(:,i)=sin(a(:,i)).*prod(cos(a(:,1:i-1)),2);
end

v(:,r-1)=prod(cos(a(:,:)),2);
v(:,r)=0;
v(:,r)=1-sum(v')';

v=v+ones(r,1)*(U1(1,:)-sum(v)/r);  %eltolja a szimplexet, hogy az U1(1,:) pont a belsejében legyen

for i=1:r
    U2=U1*inv(v);
    m=min(U2(:,i));
    v=v-m*(v-ones(r,1)*v(i,:)); %i. pontból nagyít, h. a szemközti lap érintse a ponthslmazt
    v(:,r)=0;
    v(:,r)=1-sum(v')';   
end

U2=U1*inv(v);
v;