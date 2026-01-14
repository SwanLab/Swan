sim('wing');

a = logsout.pitch;
h = logsout.plunge;
t = logsout.torque;

figure(3)
plot(a.time, a.data)
xlabel('time [s]')
ylabel('pitch [degree]')
grid on;

figure(4)
plot(h.time, h.data)
xlabel('time [s]')
ylabel('plunge [m]')
grid on;

figure(5)
plot(t.time, t.data)
xlabel('time [s]')
ylabel('torque [Nm]')
grid on;
