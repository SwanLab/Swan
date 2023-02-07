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


nVec = 1:448;
figure
plot(nVec,obj.lGtrial,nVec,obj.lGmax,nVec,obj.lG,nVec,obj.lJ,nVec,obj.l)
grid on
grid minor
xlabel('Iterations','Interpreter','latex')
ylabel('Lagrange multipliers','Interpreter','latex')
legend('$\lambda_G^{trial}$','$\lambda_G^{max}$','$\lambda_G$','$\lambda_J$','$\lambda$','Interpreter','latex')
ylim([-50 50])