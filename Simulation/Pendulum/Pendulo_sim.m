%SIMULATION PENDULUM

% Step 1: Modeling

%%Step 2: Choosing Parameters of simulation

% INITIAL CONDITIONS:
thetha_0 = 1; thetha_vel_0 = 0;
vec_InitialState = [thetha_0,thetha_vel_0];
%Constants:
Mass = 1;
L = 3;
g= 9.81;

%Simulation time, from t=0 until Tfinal (Random):
Tfinal = 10;
%Numerical integration

%%Step 3 simulation. Cicle of simulation (Euler)
%Step of simulation
h = 0.00005;  %Simulation states 0, h, 2h, 3h,...,Tfinal
%calculation of simulation
Ninstant = round(Tfinal/h)+1;

disp('Euler crono')
tic %Begining of crono
%memory reservation storing result...
time = zeros(1,Ninstant);
Result_simulation = zeros(length(vec_InitialState),Ninstant);
%Results of simulation of state "i" during instant "k" is stored in
%Result_simulation(i,k)

%X is a matrix 2 * number of points simulation
Result_simulation(:,1)=vec_InitialState; %First result, initial condition data.
time(1)=0; %Initial instant

%cycle of simulation:
for k=2:Ninstant %Virtual clock
    former_state = Result_simulation(:,k-1);
    former_time = time(k-1);
    actual_state = former_state + Model_Pendulum(former_time,former_state).*h;

    Result_simulation(:,k) = actual_state;
    time(k) = former_time +h;
end
toc

%%Results:
figure(1)
plot(time,Result_simulation), grid on
legend('Position Euler [rad]','Speed Euler [rad/s]','Location','best')
ylabel('Pos (º), Speed (rad/s)'), xlabel('time (s)')
title('Pendulum Simulation')

%%Alternative ode45
options=odeset('RelTol',1e-8);
disp('Crono ode45')
tic
[time_matlab,Result_matlab] = ode45(@Model_Pendulum,[0 Tfinal],vec_InitialState,options);
toc
%% Data procesing.
Energy = zeros(2,length(time_matlab));
Energy(1,:)=0.5.*Mass.*L.^2.*Result_matlab(:,2).^2;
Energy(2,:)=Mass.*g.*L.*(1-cos(Result_matlab(:,1)));

%%Plot alternative
hold on
plot(time_matlab,Result_matlab(:,1),'r')
plot(time_matlab,Result_matlab(:,2),'k')
hold off
legend('Position Euler [rad]','Speed Euler [rad/s]','Pos ode45','Speed ode 45','Location','best')
ylabel('Pos (º), Speed (rad/s)'), xlabel('time (s)')
title('Pendulum Simulation')


%plot(time_matlab,Energy(1,:))
%hold on
%plot(time_matlab,Energy(2,:))

%legend('Kinetic Energy','Potential Energy','Location','best')
%ylabel('Energy (W)'), xlabel('time (s)')
%title('Pendulum Simulation')