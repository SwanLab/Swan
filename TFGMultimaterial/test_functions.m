psi1 = -ones(length(p),1);
psi2 = 0*psi1;
% theta = linspace(0,2*pi,5);
% x = 0.25*cos(theta) - 0.5;
% y = 0.25*sin(theta) - 1;
% figure,plot(x,y)
hole1 = find( ((p(1,:)-0.5).^2 + (p(2,:)-1-0.125).^2)<=0.125 );
hole2 = find( ((p(1,:)-0.5).^2 + (p(2,:)-1+0.125).^2)<=0.125 );
psi2(hole1) = -1;
psi1(hole2) = 0;

chi1 = psi1>0;
chi2 = psi2>0;

mat1 = (1-chi2).*chi1;
mat2 = (chi1.*chi2);


map1 = [1, 1, 1
       1, 0, 0];
   
map2 = [1, 1, 1
       0, 0, 1];

figure(1); clf;
set(1,'WindowStyle','docked');
h1 = pdeplot(p,e,t,'xydata',-mat2,'xystyle','interp','colormap',map1,...
              'xygrid','off','colorbar','off'); axis image; axis off;
figure(2); clf;
set(2,'WindowStyle','docked');   
hold on,
         h2 = pdeplot(p,e,t,'xydata',-mat1,'xystyle','interp','colormap',map2,...
              'xygrid','off','colorbar','off'); axis image; axis off;
          hold on,
%        pdemesh(p,e,t),plot(p(1,mat2==1),p(2,mat2==1),'ks')




figure(3),clf, set(2,'WindowStyle','docked');   
% 
h1 = pdesurf(p,t,mat2);
% % colormap([0 0 1])
hold on, 
h2 = pdesurf(p,t,mat1);
hold on,
h3 = pdesurf(p,t,1-mat1-mat2); view(2)
% colormap([1 0 0])
h1.FaceColor = [1 0 0];
h2.FaceColor = [0 0 1];
h3.FaceColor = [1 1 1];







          
figure(3); clf;
set(3,'WindowStyle','docked');   
          pdeplot(p,e,t,'xydata',[-mat1,-mat2],'xystyle','flat','colormap',[map1;map2],...
              'xygrid','off','colorbar','off'); axis image; axis off;          

X = p(1,:); Y = p(2,:);
figure(4); clf;
set(5,'WindowStyle','docked');

% points = (mat1<1);
% po1 = p(:,points);
%  figure(30),clf,
%     pdeplot(po1,e,t,'xydata',mat1(points),'xystyle','interp','colormap',map1,...
%                    'xygrid','off','colorbar','on'); axis image; axis off;



    % filtering stress
    tpsi1 = pdeintrp(p,t,mat1); el1 = (tpsi1>0); to1 = t(:,el1);    
    tU1 = pdeintrp(p,t,mat1); aux1 = ones(length(tU1(el1)));    
    figure(3),clf,
    pdeplot(p,e,to1,'xydata',aux1,'xystyle','interp','colormap',map1,...
                   'xygrid','off','colorbar','on'); axis image; axis off;
               hold on, 
    tpsi2 = pdeintrp(p,t,mat2); el2 = (tpsi2>0); to2 = t(:,el2);    
    tU2 = pdeintrp(p,t,mat2); aux2 = ones(length(tU1(el2)));    
    figure(4),clf,
    pdeplot(p,e,to2,'xydata',aux2,'xystyle','interp','colormap',map2,...
                   'xygrid','off','colorbar','on'); axis image; axis off;               
               
    title('temperature'); set(1,'WindowStyle','docked');
    if(minU < maxU) 
        caxis([minU, maxU]); 
    end
