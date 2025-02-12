% Obtain filter
%a.fileName = ...;
s = FemDataContainer(a);
d.type = 'circleInclusion';
d.coord = s.mesh.coord;
d.ndim = s.mesh.ndim;
ax = 100;
ay = 1/ax;
r = 0.2185;
P = (sqrt(ax)+sqrt(ay))*(pi*r+2);
d.fracRadius = r*2;
LS = LevelSetFactory.create(d);
psi = LS.getValue();
l.value = psi;
l.type = 'LevelSet';
l.mesh = s.mesh;
LS = LevelSet(l);
f.mesh = s.mesh;
f.quadratureOrder = 'LINEAR';
f.designVariable = LS;
f.LHStype = 'AnisotropicDiffReactRobin';
chiFilter = Filter_PDE_LevelSet(f);
eh = 1:0.01:10;
ehvec = [1,3,5,10];
% eh = 10;

% Obtain perimeter
for j=1:length(eh)
    epsilon = eh(j)*s.mesh.computeMeanCellSize();
    chiFilter.updateEpsilon(epsilon);
    regularizedDensity = chiFilter.getP1fromP1(LS.value);
    regularizedDensityProjection = chiFilter.integrate_L2_function_with_shape_function(LS.value);
    Pvec(j) = sum(2/(epsilon)*((1 - regularizedDensity).*regularizedDensityProjection));
end
figure
plot(eh,Pvec,eh,P*ones(size(Pvec)),ehvec,Pvec([1,201,401,901]),'ok')
grid on
grid minor
legend('Numerical','Analytical')
xlabel('$\epsilon/h$','Interpreter','latex')
ylabel('Perimeter','Interpreter','latex')
xlim([1 10])

% Print
epsilon = eh*s.mesh.computeMeanCellSize();
chiFilter.updateEpsilon(epsilon);
regularizedDensity = chiFilter.getP1fromP1(LS.value);
ss.mesh = s.mesh;
ss.fValues = regularizedDensity;
ss.filename = 'seminarExample';
ss.order   = 'P1';
rhoe = LagrangianFunction(ss);
rhoe.print(ss);
