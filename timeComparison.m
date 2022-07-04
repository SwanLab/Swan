% e = [200,800,1250,1800,2450,3200,4050];
% SQP = [20.5,148,303.2,700];
% IPOPT = [18.7,33.57,83.93,150.9,600];
% Bisection = [22,32.3,28,24.6,30,38.3,49.41];
% AL = [50,73.4,60.2,71.6,83.2,98.4,141.92];
% Null = [25,46,88.2,93,159.7,189.8,232.3];

% figure()
% hold on
% x = length(SQP);
% plot(e(1:x),SQP,e(1:x+1),IPOPT,e,Bisection,e,AL,e,Null)
% xlabel('Number of elements')
% ylabel('Run time (s)')
% xlim([200, 4050])
% ylim([0,400])
% legend('Fmincon-SQP','Fmincon-IPOPT','Bisection','Augmented Lagrangian','Null Space')
% 
% 
% Null = [0,4209.7];

% figure()
% Sol = [300 534; 412 867; 521 902];
% x   = ['Augmented Lagrangian','Null Space', 'MMA'];
% bar(Sol)
% ylabel('Run-time (s)')
% set(gca,'xticklabel',{'Aug. Lagrangian ','Null Space','MMA'})

% load('NullSpaceBridgeC2')
% i = 1:length(c);
% subplot(1,2,1)
% plot(i,h(1,:),i,h(2,:),i,h(3,:))
% xlabel('Iterations')
% ylabel('Constraint')
% legend('g_1','g_2','g_3')
% xlim([0,i(end)])
% subplot(1,2,2)
% plot(i,d(1,:),i,d(2,:),i,d(3,:))
% xlabel('Iterations')
% ylabel('\lambda')
% legend('\lambda_1','\lambda_2','\lambda_3')
% xlim([0,i(end)])


