load("AugmentedLagrAcademic4.mat")
a.cost  = c;
a.const = h;
a.dual  = d;
a.iter  = 1:length(c);

load("NullSpaceAcademicT4.mat")
n.cost  = c;
n.const = h;
n.dual  = d;
n.iter  = 1:length(c);

load("fminconIPOPTacademic4.mat")
i.cost  = c;
i.const = h;
i.iter  = 1:length(c);

load("fminconSQPacademic4.mat")
s.cost  = c;
s.const = h;
s.iter  = 1:length(c);

subplot(1,2,1)
plot(a.iter,a.cost,'b',n.iter,n.cost,'g',i.iter,i.cost,'c',s.iter,s.cost,'r')
xlabel('Iteration')
ylabel('Cost f(x)')
legend('Augmented Lagrangian','Null Space','IPOPT','SQP')


subplot(1,2,2)
plot(a.iter,a.const(1,:),'b',n.iter,n.const(1,:),'g',i.iter,i.const(1,:),'c',s.iter,s.const(1,:),'r',...
    a.iter,a.const(2,:),'b--',n.iter,n.const(2,:),'g--',i.iter,i.const(2,:),'c--',...
    s.iter,s.const(2,:),'r--')
xlabel('Iteration')
ylabel('Constraints')
legend('Augmented Lagrangian','Null Space','IPOPT','SQP')

% figure()
% plot(a.iter,a.dual(1,:),'b',n.iter,n.dual(1,:),'r',a.iter,a.dual(2,:),'b--',n.iter,n.dual(2,:),'r--')
% xlabel('Iteration')
% ylabel('\lambda')
% xlim([1,a.iter(end)])
% legend('Augmented Lagrangian','Null Space')