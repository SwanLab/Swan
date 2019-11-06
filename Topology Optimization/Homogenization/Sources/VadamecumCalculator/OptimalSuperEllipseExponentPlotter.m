filePath = 'Topology Optimization/Vademecums/OptimalSuperEllipseExponentData.mat';
d = load(filePath);
x = d.x.thix;
y = d.x.rho;
z = d.x.q;

tri = delaunay(x,y);
ncolors = 50;
tricontour(tri,x,y,z,ncolors);

colorbar
hold on
plot(x,y,'+');
ylim([0 1])
xlabel('$\xi$','Interpreter','latex');
ylabel('\rho');
tN = '\textrm{Optimal SuperEllipse exponent}';
title(['$',tN,'$'],'interpreter','latex')
set(gca,'xtick',[0:pi/8:pi/2]) % where to set the tick marks
set(gca,'xticklabels',{'0','\pi/8','\pi/4','3\pi/8','\pi/2'})

figure

plot3(x,y,z,'+')
[xq,yq] = meshgrid(x,y);
zq = griddata(x,y,z,xq,yq);
surf(xq,yq,zq)