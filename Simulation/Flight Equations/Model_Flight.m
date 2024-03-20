function derivationofstatevector = Model_Flight(time,vec_estate,vec_constant,vec_variable)
%"normalized" form for packatges of the simulation:
%Defining "state vector"
%vec_estate=[theta, theta_vel]

x = vec_estate(1);
hight = vec_estate(2);
v = vec_estate(3);
gamma = vec_estate(4);
Weight = vec_estate(5);

%Constant parameters of the system:


cT = vec_constant(1);
S = vec_constant(2);
g = vec_constant(3); 
T_0 = vec_constant(4);
R = vec_constant(5);

%Vec_variable:

Pow_set = vec_variable(1);
alpha = vec_variable(2);


%Entry: 

if hight<11000
    dT = -0.0065;
    T = 288.15+dT.*(hight-0);
    rho = 1.225.*(T./288.15).^(-1-g./(dT.*R));
    
elseif hight<=20000
    T = 216.65;
    rho = 0.3639.*exp((-g./(R.*T_0)).*(hight-11000));
    
else
    T = 216.65;         %Just in case h>= 20.000 m
    rho = 0.0880;
end


a = 331.*sqrt(T./273);  %Speed of sound.
Mach = v./a;


Cd_0 = 0.04;    %Modificar dependiendo de mach
mu = 0.1;
C_L_alpha = 4;



Cd = Cd_0 + mu.*C_L_alpha.*alpha.^2;
Cl = C_L_alpha.*alpha;

Drag = 0.5.*Cd.*S.*rho.*v.^2;
Lift = 0.5.*Cl.*S.*rho.*v.^2;
%Thrust = (rho./1.225)
Thrust = (rho./1.225).*22.*10e3.*Pow_set;



%Elemental equations:
dx = v.*cos(gamma);
dh = v.*sin(gamma);
v_dt = g.*((Thrust-Drag)./Weight-sin(gamma));
d_gamma = (g./v).*(Lift./Weight-cos(gamma));
d_Weight = -cT.*Thrust;

derivationofstatevector=...
    [dx;...    
     dh;...
     v_dt;...    
     d_gamma;...
     d_Weight];

end

