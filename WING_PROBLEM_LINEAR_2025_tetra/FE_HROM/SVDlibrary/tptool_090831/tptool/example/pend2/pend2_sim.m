sim('pend2');

a = logsout.angles;
f = logsout.force;

figure(5)
plot(a.time, a.data*180/pi)
xlabel('time')
ylabel('angle [degree]')
grid on;

figure(6)
plot(f.time, f.data)
xlabel('time')
ylabel('force [N]')
grid on;
