% TESTING CHARACTERISTIC FUNCTION
%% Get mesh
clear; close all;

file = 'RVE_Square_Triangle';
% file = 'BridgeCool_Quadrilateral_Bilinear_Structured';
a.fileName = file;
m = FemDataContainer(a);

%% Create characteristic function
s.mesh    = m.mesh;
s.fxy     = @(x,y) (x-0.5)^2+(y-0.5)^2-0.3^2;
% s.fxy     = @(x,y) (x-1.5)^2+(y-1)^2-0.6^2;
circleFun = CharacteristicFunction(s);

%% L2 Projectors
% CharFun to P1 Function
x.mesh   = m.mesh;
x.connec = m.mesh.connec;
projP1 = Projector_toP1(x);
resCharFunInP1 = projP1.project(circleFun);
resCharFunInP1.plot(m.mesh);
% view(0,90)
title('L2p1')

% CharFun to P0 Function
projP0 = Projector_toP0(x);
resCharFunInP0 = projP0.project(circleFun);
resCharFunInP0.plot(m.mesh);
title('L2p0')


% CharFun to P1 Discontinuous Function
projP1D = Projector_toP1Discontinuous(x);
resCharFunInP1D = projP1D.project(circleFun);
resCharFunInP1D.plot(m.mesh);
title('L2p1d')

%% H1 Projectors

% CharFun to P1 Function in H1 domain
xx = x;
xx.projectorType = 'toH1P1';
xx.epsilon = m.mesh.computeMeanCellSize();
projH1P1 = Projector.create(xx);
resCharFunInH1P1 = projH1P1.project(circleFun);
resCharFunInH1P1.plot(m.mesh);
% view(0,90)
title('H1p1')

% CharFun to P1 Discontinuous Function in H1 domain
xx.projectorType = 'toH1P1Disc';
projH1P1D = Projector.create(xx);
resCharFunInH1P1D = projH1P1D.project(circleFun);
resCharFunInH1P1D.plot(m.mesh);
title('H1p1d')

%% Filters
% Heaviside
coor = m.mesh.coord;
fxy = zeros(size(coor,1),1);
for i = 1:size(coor,1)
    fxy(i) = s.fxy(coor(i,1),coor(i,2));
end

% Design variable
fileName = 'test_FilterVSProjectors';
q = Settings(fileName);
q.warningHoleBC = false;
q.printIncrementalIter = false;
q.printChangingFilter = false;
q.printing = false;
translator = SettingsTranslator();
translator.translate(q);
fileName = translator.fileName;
settings  = SettingsTopOptProblem(fileName);
settings.designVarSettings.initialCase = 'given';
settings.designVarSettings.creatorSettings.value = fxy;
q = settings.designVarSettings;
q.mesh = Mesh_Total(m.mesh);
q.scalarProductSettings.epsilon = m.mesh.computeMeanCellSize();
q.scalarProductSettings.mesh = m.mesh;
ls = LevelSet(q);

% Filter
y = SettingsFilter();
y.filterType = 'PDE';
y.quadratureOrder = 'LINEAR';
y.designVarType = 'LevelSet';
y.femSettings.scale = 'MACRO';
y.femSettings.mesh = m.mesh;
y.mesh = m.mesh;
y.designVariable = ls;
filter = Filter.create(y);
filteredVals = filter.getP1fromP1(fxy);
z.connec  = m.mesh.connec;
z.type    = m.mesh.type;
z.fValues = filteredVals;
filterToP1 = P1Function(z);
filterToP1.plot(m.mesh);
title('PDE Filter');

% 1D - PLOT in function of theta; value(nodes)=f(x)
%    where nodes: 0-eps <= y-mx <= 0+eps
X = m.mesh.coord(:,1)-0.5;
Y = m.mesh.coord(:,2)-0.5;
theta = 0;
slope = tan(theta);
eps = m.mesh.computeMeanCellSize();
fun = Y-slope*X;
k1 = fun >= 0-eps;
k2 = fun <= 0+eps;
nodesk = find(k1==k2);
l = sqrt(X(nodesk).^2+Y(nodesk).^2);
figure
plot(l,filterToP1.fValues(nodesk),'.')
hold on
plot(l,resCharFunInH1P1.fValues(nodesk),'.')
hold on
plot(l,resCharFunInP1.fValues(nodesk),'.')
hold off
grid on
grid minor
xlabel('Length','Interpreter','latex')
ylabel('CharFun','Interpreter','latex')
legend('PDE Filter','P1(H1) Projection','P1(L2) Projection','Interpreter','latex')
ylim([-0.1 1.1])

%% Perimeter evolution vs epsilon
chi = circleFun;

xx = x;
xx.projectorType = 'toH1P1';
epsOpt = m.mesh.computeMeanCellSize();
eps = 0.5*epsOpt:0.25*epsOpt:100*epsOpt;
P = zeros(size(eps));

for k=1:length(eps)
    xx.epsilon = eps(k);
    projH1P1 = Projector.create(xx);
    rhoe = projH1P1.project(circleFun);
    P(k) = computePerimeterIntegral(chi,rhoe,eps(k),m);
end

Ptheo = 2*pi*0.3*ones(size(P));

figure
plot(eps/epsOpt,Ptheo,eps/epsOpt,P)
grid on
grid minor
legend('Theoretical value','Numerical value','Interpreter','Latex')
xlabel('$\epsilon/h$','Interpreter','Latex')
ylabel('Perimeter','Interpreter','Latex')

%% Functions

function P = computePerimeterIntegral(chi,rhoe,eps,m)

order = 'QUADRATIC';
q = Quadrature.set(m.mesh.type);
q.computeQuadrature(order);
xV = q.posgp;
dV = m.mesh.computeDvolume(q);
nGaus = q.ngaus;
nElem = m.mesh.nelem;

fChi = chi.evaluate(xV);
fRhoe = rhoe.evaluate(xV);
P = 0;

for i = 1:nElem
    for igaus = 1:nGaus
        dVg = dV(igaus, i);
        fC = fChi(1,igaus,i);
        fR = fRhoe(1,igaus,i);
        P = P + (1-fR)*fC*dVg;
    end
end

P = P/(2*eps);

end