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