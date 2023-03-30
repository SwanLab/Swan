m=obj.designVariable.getUnfittedMesh
s.connec = m.backgroundMesh.connec;
s.type = m.backgroundMesh.type;
s.fValues = DmF;
P1fun = P1Function(s);
P1fun.plot(m.backgroundMesh)



m=obj.designVariable.getUnfittedMesh
figure
m.plot


r = 1.1;
n=1:k;
tau=1./(r.^(n-1));
figure
plot(tau,obj.meritEv,tau,obj.meritOldEv,'red')
grid on
grid minor
xlabel('$\tau$','Interpreter','latex')
ylabel('Merit function','Interpreter','latex')
legend('Current','Old','Interpreter','latex')


nVec = 1:length(obj.lG);
figure
plot(nVec,obj.lG,nVec,obj.lJ)
hold on
plot(nVec,obj.l,'Color',"#EDB120")
grid on
grid minor
xlabel('Iterations','Interpreter','latex')
ylabel('Lagrange multipliers','Interpreter','latex')
legend('$\lambda_G$','$\lambda_J$','$\lambda$','Interpreter','latex')
ylim([-52 4])

figure
semilogy(nVec,obj.aJvec,nVec,obj.aGvec)
grid on
grid minor
xlabel('Iterations','Interpreter','latex')
ylabel('Nullspace coefficients','Interpreter','latex')
legend('$a_J$','$a_G$','Interpreter','latex')
ylim([-0.1 1])









% Obtain filter
a.fileName = 'SquareForAniTests';
s = FemDataContainer(a);
d.type = 'circleInclusion';
d.coord = s.mesh.coord;
d.ndim = s.mesh.ndim;
r = 0.2185;
P = 2*pi*r+4;
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
f.LHStype = 'DiffReactRobin';
chiFilter = Filter_PDE_LevelSet(f);
eh = 1:0.01:10;
eh = 10;

% Obtain perimeter
for j=1:length(eh)
    epsilon = eh(j)*s.mesh.computeMeanCellSize();
    chiFilter.updateEpsilon(epsilon);
    regularizedDensity = chiFilter.getP1fromP1(LS.value);
    regularizedDensityProjection = chiFilter.integrate_L2_function_with_shape_function(LS.value);
    Pvec(j) = sum(2/(epsilon)*((1 - regularizedDensity).*regularizedDensityProjection));
end
figure
plot(eh,Pvec,eh,P*ones(size(Pvec)))
grid on
grid minor
legend('Numerical','Analytical')
xlabel('$\epsilon/h$','Interpreter','latex')
ylabel('Perimeter','Interpreter','latex')

% Print
epsilon = eh*s.mesh.computeMeanCellSize();
chiFilter.updateEpsilon(epsilon);
regularizedDensity = chiFilter.getP1fromP1(LS.value);
ss.mesh = s.mesh;
ss.fValues = regularizedDensity;
ss.filename = 'seminarExample';
rhoe = P1Function(ss);
rhoe.print(ss);