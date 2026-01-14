%k�zel-NO dekompoz�ci� r-1 dimenzi�ra
%
% [U2,fi2]=CNO(U1,h,hh,n1,n2,n3));
%
%az U1 n-szer r-es SN tipusu m�trix,n>=r, r>=2
%U2 (m�rete nxr) �s fi2 (m�rete rxr) is SN, U2*fi2=U1
%�s U2 k�zel-NO, azaz oszlopainak legnagyobb eleme majdnem 1.
%
%Haszn�lja a closeness.m �s a polarorto.m f�jlt.
%

% %PARAMETEREK
% h: kb. 0.5 �s 5 k�z�tt �rdemes megv�lasztani. Kis h eset�n n�h�ny 
% oszlopmaximum nagyon k�zel lesz 1-hez, nagy h eset�n az oszlopmaximumok
% kb. egyenl? t�vra lesznek 1-t?l. tul kis vagy tul nagy h1 eset�n 
% numerikus hibak lephetnek fel
%
%hh: hh s�llyal veszi az optimum keres�s�ekor figyelembe azt, hogy a 
%f�ggv�nyg�rb�k v�gei (bal �s jobb) k�z�tt van e 1-hez k�zeli 
%hh-t 0 �s kb. r k�z�tt c�lszer? megv�lasztani.  
% 
% (n1, n2, n3): eg�sz sz�mok min�l nagyobbak, ann�l pontosabb
% de lassabb a fut�s. Erdemes tobbszor futtatni �s ha nagyon elt�r? 
% eredm�nyeket ad, e h�rom param�ter valamelyik�t, elsosorban n1-et 
% novelni. Ajanlott ertekek r=4 eseten kb. (10-30,0-5,1-5).

%A f�ggv�ny r=3 eset�n n�h�ny m�sodpercig, r=4 eset�n kb 
%1 percig fut. A vegen kiir egy "NOtol valo tavolsag"-ot 
%jellemzo meroszamot.
%


function [U2,fi2] = CNO(U1,h,hh,nveletlen,nzavar,nlokalis);

%global U1

% nveletlen=10;
% nzavar=4;

% path(path,'c:\peti\baranyip');
%path(path,'d:\closet\ndimenzio');


[n r] = size(U1);
if r>n
	error('Ha a bemeno matrix nxr-es, r>n nem lehet');
end
%if norm(1-sum(U1'))>10^(-6)
%    error('A bemeno matrix nem SN')
%end

if r==2
	[fi2 i] = min(U1); %a fi nem kell semmire
	fi2 = [U1(i(1),:); U1(i(2),:)];
	U2 = U1/fi2;
else
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      ponthalmaz inercia
	center = ones(n,1)*sum(U1(:,1:r-1)/n);
	U11 = U1(:,1:r-1)-center;
	inercia = U11'*U11;
	%%%%%%%%%%%%%%%%% transzform�ci�, hogy egys�ginerci�ja legyen
	[sv,se] = eig(inercia);
	U111 = U11*sv*se^(-.5);
	%plot3(U111(:,1),U111(:,2),U111(:,3),'r')
	%plot(U111(:,1),U111(:,2))
	
	U1 = U111;                                %transzform�lt m�trixszal folytatjuk!!!
	U1(:,r) = 1-sum(U1')';
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     idealis felbontas numerikus uton 
	
	ertekmin = 1000;
	%%%%%%%%%keres veletlen helyeken
	rand('state',sum(100*clock));
	for i = 1:nveletlen;
		i;
		veletlenhely=(rand(r,r-2)-.5)*pi; %r-1D-s polarkoordinatak (r-2 db.) a -pi/2...pi/2 intervallumban
		veletlenhely(:,r-2)= veletlenhely(:,r-2)*2; %az utolso polarkoordinata -pi...pi
		[polar ertek] = fminsearch(@closeness,veletlenhely,[],U1,h,hh);
		
		
		if ertek<ertekmin   %ha az eddigi legjobb megoldast talalta...
			ertekmin=ertek;
			polarmin=polar;
			disp(ertek);
			%warning('siker- veletlen probalkozas')
		end
	end
	
	%%%%%%%%%%a legjobb megoldast lokalisan javitja
	szamlalo=0;
	while szamlalo < nlokalis
		[polar,ertek] = fminsearch(@closeness,polarmin,[],U1,h,hh);
		if ertek < ertekmin-.001
			szamlalo=0;
			ertekmin=ertek;
			polarmin=polar;
			%warning('siker- direkt javitas')
		else
			szamlalo=szamlalo+1;
		end
	end
	disp(ertek);
	
	%%%%%%%%%%a legjobb megoldas kis zavarasa
	for i=1:nzavar
		i;  
		zavar=(rand(r,r-2)-.5)*pi;
		zavar(:,r-2)=zavar(:,r-2)*2;
		[polar,ertek] = fminsearch(@closeness,polarmin+.1*rand(1)*zavar,[],U1,h,hh);
		%%%%%%%%%a megoldast lokalisan javitja
		szamlalo=0;
		while szamlalo<nlokalis
			[polarseg,ertekseg] = fminsearch(@closeness,polar,[],U1,h,hh);
			if ertekseg < ertek-.001
				szamlalo=0;
				%warning('siker- direkt javitas')
				ertek=ertekseg;
				polar=polarseg;
			else
				szamlalo=szamlalo+1;
			end
		end
		if ertek<ertekmin
			disp(ertek);
			ertekmin=ertek;
			polarmin=polar;
			%warning('siker-utozavaras');
		end
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%az eredm�ny egys�ges�t�se
	if r>3
		for i=1:r
			for j=1:r-3
				if mod(polarmin(i,j)+pi/2,2*pi)>pi
					polarmin(i,j+1)=polarmin(i,j+1)+pi;
					polarmin(i,j)=pi-polarmin(i,j);
				end
			end
		end
		polarmin(:,r-3)=mod(polarmin(:,r-3),2*pi);
	end
	polarmin(:,1:r-2)=mod(polarmin(:,1:r-2)+pi/2,2*pi)-pi/2;
	
	polarmin=sortrows(polarmin);
	
	disp('Distance from NO')
	disp(ertekmin)
	
	[U2 fi2]=polarorto(U1,polarmin);
	
	
	
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%az eredm�ny �br�zol�sa r=3,4 eset�n
%    if r==4
%        figure(9);
%        hold on;
%        plot3(U1(:,1),U1(:,2),U1(:,3),'r')
%        for i=1:4
%            for j=1:i-1
%                plot3([fi2(i,1) fi2(j,1)],[fi2(i,2) fi2(j,2)],[fi2(i,3) fi2(j,3)]);
%            end
%        end
%    end
	
%    if r==3
%        figure(9);
%        hold on;
%        plot(U1(:,1),U1(:,2),'r')
%        for i=1:3
%            for j=1:i-1
%                plot([fi2(i,1) fi2(j,1)],[fi2(i,2) fi2(j,2)]);
%            end
%        end
%    end
	
	%visszatranszform�l�s
	fi2=fi2(:,1:r-1)*se^.5*sv'+center(1:r,:);
	fi2(:,r)=0;
	fi2(:,r)=1-sum(fi2')';
end

