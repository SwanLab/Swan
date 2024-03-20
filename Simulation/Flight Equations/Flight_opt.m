%Bisection method:

%%Input argument:
%v_0_max = 1000; %Max initial speed [m/s]
%v_0_min = 1; %Min initial speed [m/s]
%Tolerance = 10e0;%Tolerance of the numerical result.
%Max_Iter = 80; %Max numer of iteration.
h = 0.001; %Step simulation.
%Constants
L = 2000;

%theta_0 = pi./4;
rho = 1.225; %Density kg/m3

cT = 0.00001;
S = 34;
g = 9.81;
T_0 = 288.15;
R = 287;
vec_constant = [cT,S,g,T_0,R];
%Boundary conditions:
x_0 = 0;
hight_0 = 0;
v_0 = 100; %m/s
gamma_0 = 0.1; %10º;
Weight = 6000.*g; %Kg

%Variable set
Pow_set =0.7;
Alpha = 0.02;
vec_variable = [Pow_set, Alpha];



%%Loop iteration:
%N = 1;
%v_a = v_0_max;
%v_b = v_0_min;
disp('Crono:')
tic

%while N<= Max_Iter  %Limiting iteractions to prevent infinte loop.
    [Result_sim_C,Hit_time_C,Time] = Flight_sim(vec_constant,x_0,hight_0,v_0,gamma_0,Weight,vec_variable,h);
    %diff_C = Result_sim_C(1,int32(Hit_time_C/h+1))-vec_constant(1);
   % if diff_C >= -Tolerance/2 && diff_C <= Tolerance/2
        Final_Result = Result_sim_C;
        Final_Hit_time = Hit_time_C;
        %break
    %end
   % if sign(diff_C)==sign(diff_A)
       % v_a = v_c;
    %elseif sign(diff_C)==sign(diff_B)
        %v_b = v_c;
    %end
    %disp(N)
    %N = N+1;
%end

toc
k_final = int32(Final_Hit_time/h+1);

figure(1);
plot(Final_Result(1,1:k_final),Final_Result(2,1:k_final))
xlabel('X [m]'); ylabel('Y[m]');
title('Launch Simulation, y(x)')

figure(2)
plot(Time(1:k_final),Final_Result(1,1:k_final))
xlabel('t [s]'); ylabel('X [m]');
title('Launch Simulation, x(t)')

figure(3)
plot(Time(1:k_final),Final_Result(2,1:k_final))
xlabel('t [s]'); ylabel('Y [m]');
title('Launch Simulation, y(t)')

figure(4)
plot(Time(1:k_final),Final_Result(3,1:k_final))
xlabel('t [s]'); ylabel('Speed [m/s]');
title('Launch Simulation, dt(t)')

figure(5)
plot(Time(1:k_final),Final_Result(4,1:k_final))
xlabel('t [s]'); ylabel('Gamma [rad]');
title('Launch Simulation, Gamma(t)')
