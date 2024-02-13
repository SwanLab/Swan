function derivationofstatevector = Model_Lanz(time,vec_estate,vec_constant)
%"normalized" form for packatges of the simulation:
%Defining "state vector"
%vec_estate=[theta, theta_vel]

x = vec_estate(1);
y = vec_estate(2);
dx = vec_estate(3);
dy = vec_estate(4);
%Constant parameters of the system:


g = vec_constant(2);

%Entry: 

%Elemental equations:


derivationofstatevector=...
    [dx;...    
     dy;...
     0;...    
     -g];

end

