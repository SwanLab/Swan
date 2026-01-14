sim('pend');

a = logsout.angle;
f = logsout.force;

figure(3)
plot(a.time, a.data*180/pi)
xlabel('time')
ylabel('angle [degree]')
grid on;

figure(4)
plot(f.time, f.data)
xlabel('time')
ylabel('force [N]')
grid on;
