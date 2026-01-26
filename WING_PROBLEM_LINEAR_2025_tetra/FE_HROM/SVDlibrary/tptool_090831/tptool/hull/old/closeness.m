%U1nxr-es mátrix lehet,SN
%megadja az adott alakú (polarv) INO burkolószimplex 
%távolságát az NO feltételtol

%norma=closeness(polarv,U1,h)

%polarv a tetraéder körülírt köre közepéböl a 
%csúcsokba mutató vektorok iránya r-1 ndimenzios polárkoordinátákkal
%norma pedig 1-az oszlopmaximumok h-normája
%ill ha hh==1, még hozzáveszi az els? és az utols? sor maximumait 

function norma=closeness(polarv,U1,h,hh)

[n,r]=size(U1);

if size(polarv)~=[r r-2]
	error('ha U1 nxr-es, akkor polarv rx(r-2)-es kell legyen')
end
if n<r
	error('ha U1 mérete nxr, akkor n>=r kell teljesüljön')
end

a=polarv;

for i=1:r-2	
	v(:,i)=sin(a(:,i)).*prod(cos(a(:,1:i-1)),2);
end
v(:,r-1)=prod(cos(a(:,:)),2);
v(:,r)=1-sum(v')';

v=v+ones(r,1)*(U1(1,:)-sum(v)/r);  %eltolja a szimplexet, hogy az U1(1,:) pont a belsejében legyen

for i=1:r
	U2=U1/v;
	m=min(U2(:,i));
	v=v-m*(v-ones(r,1)*v(i,:)); %i. pontból nagyít, h. a szemközti lap érintse a ponthslmazt
	v(:,r)=0;
	v(:,r)=1-sum(v')';
end
U2=U1/v;
v;
%norma=norm(1-max(U2),h);

norma=norm(log(max(U2)),h);
vegek=log(sort(max([U2(1,:);U2(n,:)])));
norma=norma+hh*norm(vegek(r-1:r),h);

%figure(1);
%hold on;
%plot3(U1(:,1),U1(:,2),U1(:,3))
%plot3(v(:,1),v(:,2),v(:,3))
