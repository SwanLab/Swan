function [b,p,e,t]=raccommode(b1,p1,e1,t1,b2,p2,e2,t2)

% cherche les noeuds communs
np1=size(p1,2); np2=size(p2,2);
p1x=p1(1,:); p1y=p1(2,:);
p2x=p2(1,:); p2y=p2(2,:);
%figure(1); clf; pdemesh(p1,e1,t1);
%figure(2); clf; pdemesh(p2,e2,t2);
pp1x=ones(np2,1)*p1x; pp1y=ones(np2,1)*p1y;
pp2x=ones(np1,1)*p2x; pp2y=ones(np1,1)*p2y;
ddx=pp2x'-pp1x; ddy=pp2y'-pp1y;
dd=sqrt(ddx.^2+ddy.^2);
[pcom2,pcom1]=find(dd==0);
npcom=length(pcom1);

% cherche les bords en commun
[pecom11,ecom11,rien]=intersect(e1(1,:),pcom1);
[pecom12,ecom12,rien]=intersect(e1(2,:),pcom1);
ecom1=intersect(ecom11,ecom12);
[pecom21,ecom21,rien]=intersect(e2(1,:),pcom2);
[pecom22,ecom22,rien]=intersect(e2(2,:),pcom2);
ecom2=intersect(ecom21,ecom22);
necom=length(ecom1);

% cherche colonne de b en commun
bcom1=e1(5,ecom1);
bcom2=e2(5,ecom2);

% faire b
b=[b1,b2];
nb1=size(b1,2); nb2=size(b2,2);
%[b,ib1,ib2]=union(b1,b2,'rows');

% faire p
pncom2=setdiff([1:np2],pcom2);
p=[p1,p2(:,pncom2)];
pnum2=[np1+1:np1+np2];
for i=1:npcom
    pnum2(pcom2(i))=pcom1(i);
    pnum2(pcom2(i)+1:np2)=pnum2(pcom2(i)+1:np2)-1;
end;

% faire e
e2purge=e2;
e2purge(:,ecom2)=[];
e2renum=e2purge;
e2renum(1,:)=pnum2(e2purge(1,:)); e2renum(2,:)=pnum2(e2purge(2,:)); e2renum(5,:)=e2purge(5,:)+nb1;
e=[e1,e2renum];

%faire t
t2renum=t2;
t2renum(1,:)=pnum2(t2(1,:));
t2renum(2,:)=pnum2(t2(2,:));
t2renum(3,:)=pnum2(t2(3,:));
t=[t1,t2renum];



