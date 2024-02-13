function derivationofstatevector = Model_Pendulum(time,vec_estate)
%"normalized" form for packatges of the simulation:
%Defining "state vector"
%vec_estate=[theta, theta_vel]

thetha = vec_estate(1);
thehta_vel = vec_estate(2);

%Constant parameters of the system:

g = 9.81;   %Gravity m/s2
Mass = 1;   %Mass kg
L = 3;    %Lonitud [m]

%Entry: 

%Elemental equations:
%T = Mass.*L.*thehta_vel.^2+Mass.*g.*cosd(thetha);

derivationofstatevector=...
    [thehta_vel;...    
    -(g./L).*sin(thetha)];
%sqrt((T-Mass.*g.*sind(thetha))./Mass.*L)
end

