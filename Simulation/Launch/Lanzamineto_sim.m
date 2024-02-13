function [Result_simulation,Hit_time,time] = Lanzamineto_sim(vec_constant,v_0,x_0,y_0,h)
%Constants:
L = vec_constant(1);
g = vec_constant(2);
theta_0 = vec_constant(3);
v_0_x = v_0.*cos(theta_0);
v_0_y = v_0.*sin(theta_0);
vec_InitialState = [x_0,y_0,v_0_x,v_0_y];
%Simulation time, from t=0 until Tfinal (Random):
Tfinal = 1000;   
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
    actual_state = former_state + Model_Lanz(former_time,former_state,vec_constant).*h;
    Result_simulation(:,k) = actual_state;
    time(k) = former_time +h;
    if Result_simulation(2,k)<=0
        Hit_time = time(k);
        break
    end
end

end

