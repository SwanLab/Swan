sim('tora');

a = logsout.angle;
da = logsout.dangle;
t = logsout.torque;

figure(3)
plot(a.time, a.data)
xlabel('time [s]')
ylabel('angle [degree]')
grid on;

figure(4)
plot(da.time, da.data)
xlabel('time [s]')
ylabel('dangle/dt [degree/s]')
grid on;

figure(5)
plot(t.time, t.data)
xlabel('time [s]')
ylabel('torque [Nm]')
grid on;
