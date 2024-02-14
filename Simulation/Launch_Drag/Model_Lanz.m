function derivationofstatevector = Model_Lanz(time,vec_estate,vec_constant)
%"normalized" form for packatges of the simulation:
%Defining "state vector"
%vec_estate=[theta, theta_vel]

x = vec_estate(1);
y = vec_estate(2);
dt = vec_estate(3);
gamma = vec_estate(4);
%Constant parameters of the system:


g = vec_constant(2);
K = vec_constant(4);
m = vec_constant(5); 

%Entry: 

%Elemental equations:
dx = dt.*cos(gamma);
dy = dt.*sin(gamma);
dt_dt = -g.*sin(gamma)-(K./m).*dt.^2;
d_gamma = -(g./dt).*cos(gamma);

derivationofstatevector=...
    [dx;...    
     dy;...
     dt_dt;...    
     d_gamma];

end

