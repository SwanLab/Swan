function Data_representation(y,t,phisical_cons, output)
% First plot trajectory:
figure(1)
plot(y(:,1),y(:,2))
title Trajectory
xlabel [m]
ylabel [m]
axis equal

% Second plot Speed over time: 
figure(2)
plot(t(:,1),y(:,3))
title Speed
xlabel [s]
ylabel [m/s]

% Third plot Angle over time:
figure(3)
plot(t(:,1),y(:,4))
title Gamma
xlabel [s]
ylabel [rad]

gamma0 = y(1,4).*(180./pi);
X = ['Initial Angle: ' num2str(gamma0) 'º'];
disp(X)
X = ['Time: ' num2str(t(end,1)) ' s'];
disp(X)
disp(output)
end

