clear;
clc;
close all;

% Mesh
x1       = linspace(0,2,200);
x2       = linspace(0,2,200);
[xv,yv]  = meshgrid(x1,x2);
[F,V]    = mesh2tri(xv,yv,zeros(size(xv)),'x');
s.coord  = V(:,1:2);
s.connec = F;
mesh     = Mesh(s);

wVec = 0:5:90;
beta = [45,60,89,89.9];

for i = 1:length(beta)
    for j = 1:length(wVec)
        Pval(i,j) = obtainAPValue(mesh,wVec(j),beta(i));
    end
end


wVecIdeal = 0:0.01:90;

figure
plot(wVec,Pval(1,:)/Pval(1,1),wVec,Pval(2,:)/Pval(2,1),wVec,Pval(3,:)/Pval(3,1),wVecIdeal,[ones(length(wVecIdeal(wVecIdeal<45)),1);zeros(length(wVecIdeal(wVecIdeal>=45)),1)],'-k');
grid minor
xlabel('Bar angle ($^\circ$)','Interpreter','latex')
ylabel('P/P($0^\circ$)','Interpreter','latex')
legend('$\beta=45^\circ$ (isotropic)','$\beta=60^\circ$','$\beta=89^\circ$','Ideal','Interpreter','latex')
ylim([-0.1 1.1])





function Pval = obtainAPValue(mesh,w,beta)
% Bar
gPar.type          = 'RectangleRotated';
gPar.xSide         = 1;
gPar.ySide         = 0.1;
gPar.xCoorCenter   = 1;
gPar.yCoorCenter   = 1;
gPar.omegaDeg      = w;
g                  = GeometricalFunction(gPar);
lsFun              = g.computeLevelSetFunction(mesh);

ss.filterType = 'LUMP';
ss.mesh       = mesh;
ss.trial      = P1Function.create(mesh,1);
filter        = Filter.create(ss);
lsFun         = filter.compute(lsFun,'QUADRATICMASS');

sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh();
uMesh              = UnfittedMesh(sUm);
uMesh.compute(lsFun.fValues);
sss.uMesh = uMesh;
charFun = CharacteristicFunction.create(sss);

% PDE Filter
ss.filterType   = 'PDE';
ss.mesh         = mesh;
ss.boundaryType = 'Neumann';
ss.metric       = 'Anisotropy';
ss.trial        = P1Function.create(mesh,1);
ss.CAnisotropic = [tand(beta), 0; 0, 1/tand(beta)];
ss.aniAlphaDeg  = 90;
filter          = Filter.create(ss);

% Create perimeter functional
shF.mesh = mesh;
shF.filter = filter;
shF.epsilon = mesh.computeMeanCellSize();
P        = ShFunc_Perimeter(shF);
Pval     = P.computeFunction(charFun);
end