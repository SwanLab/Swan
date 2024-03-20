function [Result_simulation,Hit_time,time] = Flight_sim(vec_constant,x_0,hight_0,v_0,gamma_0,Weight,vec_variable,h)
%Constants:

L = vec_constant(1);
g = vec_constant(2);
theta_0 = vec_constant(3);


vec_InitialState = [x_0,hight_0,v_0,gamma_0,Weight];
%Simulation time, from t=0 until Tfinal (Random):
Tfinal = 3000;   
%Numerical integration

%%Step 3 simulation. Cicle of simulation (Euler)
%Step of simulation
%Simulation states 0, h, 2h, 3h,...,Tfinal
%calculation of simulation
Ninstant = round(Tfinal/h)+1;
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
    actual_state = former_state + Model_Flight(former_time,former_state,vec_constant,vec_variable).*h;
    Result_simulation(:,k) = actual_state;
    time(k) = former_time +h;
    if Result_simulation(2,k)>=2000
        Hit_time = time(k);
        break
    end
end

end

