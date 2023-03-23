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



a.fileName = 'SquareForAniTests';
s = FemDataContainer(a);
d.type = 'squareInclusion';
d.coord = s.mesh.coord;
d.ndim = s.mesh.ndim;
d.widthSquare = sqrt(0.15);
LS = LevelSetFactory.create(d);
psi = LS.getValue();
chi = 1-heaviside(psi);
ss.mesh = s.mesh;
ss.fValues = chi;
Square = P1Function(ss);
ss.filename = 'seminarExample';
% Square.print(ss);
SquareRelIso = Square.project('H1P1');
SquareRelIso.print(ss)