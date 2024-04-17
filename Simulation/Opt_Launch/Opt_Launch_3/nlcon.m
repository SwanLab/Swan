function [c, ceq,gradc,gradceq] = nlcon(x,num_cons,phisical_cons)

%Constants:
g = phisical_cons(1,1);
v_0 =phisical_cons(1,2);
x_1_0 =phisical_cons(1,3);
x_2_0 =phisical_cons(1,4);
t_0 =phisical_cons(1,5);
h = num_cons;


%inequality constrains
c = [];
% Initialize ceq as an array
    ceq = zeros(4.*h+1, 1);
%Constrains:
for i= 1:(h-1)

    %Variables before
    x_1_ant =  x(4.*(i-1)+1); 
    x_2_ant =  x(4.*(i-1)+2);
    v_ant   =  x(4.*(i-1)+3);
    gamma_ant =x(4.*(i-1)+4);
    
    %Variables after
    x_1_new =  x(4.*i+1); 
    x_2_new =  x(4.*i+2);
    v_new   =  x(4.*i+3);
    gamma_new =x(4.*i+4);

    %Diferential of variable:
    dx_1 = v_ant.*cos(gamma_ant);
    dx_2 = v_ant.*sin(gamma_ant);
    d_v = -g.*sin(gamma_ant);
    d_gamma = -(g./v_ant).*cos(gamma_ant);

    %Increment of time:
    t_f = x(4.*h+1);

    Inc_time = ((t_f-t_0)./h);


    ceq(4.*i+1,1) = x_1_ant + Inc_time.*dx_1-x_1_new;
    ceq(4.*i+2,1) = x_2_ant + Inc_time.*dx_2-x_2_new;
    ceq(4.*i+3,1) = v_ant   + Inc_time.*d_v-v_new;
    ceq(4.*i+4,1) = gamma_ant +Inc_time.*d_gamma-gamma_new;

end

%First set
    gamma_0 = x(4.*h+2);

    %Time variables
    t_f = x(4.*h+1);

    Inc_time = ((t_f-t_0)/h);
    
    %Diferentials 
    dx_1 = v_0.*cos(gamma_0);
    dx_2 = v_0.*sin(gamma_0);
    d_v = -g.*sin(gamma_0);
    d_gamma = -(g./v_0).*cos(gamma_0);

    %New variables
    x_1 = x(1);
    x_2 = x(2);
    v_1 = x(3);
    gamma_1 = x(4);


ceq(1,1) = x_1_0 + Inc_time.*dx_1-x_1;
ceq(2,1) = x_2_0 + Inc_time.*dx_2-x_2;
ceq(3,1) = v_0 + Inc_time.*d_v-v_1;
ceq(4,1) = gamma_0 + Inc_time.*d_gamma-gamma_1;
%Second set
x_2_end =  x(4.*(h-1)+2);

ceq(4.*h+1,1) = x_2_end;

%Gradient calculation:
gradc = [];
gradceq = zeros(4.*h+1,4.*h+2); % For future code constrains in colums variables in rows
%Secction 1:
for i = 1:(h-1)
    
    %Variables
    x_1_ant = x(4.*(i-1)+1);
    x_2_ant = x(4.*(i-1)+2);
    v_ant = x(4.*(i-1)+3);
    gamma_ant = x(4.*(i-1)+4);

    %Time variables:
    t_f = x(4.*h+1);

    Inc_time = ((t_f-t_0)/h);
    
    %C1
    %Diferentials
    %dgamma_gamma = (1-(g./v_ant).*sin(gamma_ant))^(-1);
    dv = cos(gamma_ant);
    dgamma = -v_ant.*sin(gamma_ant);
    dt_f = (1./h).*v_ant.*cos(gamma_ant);
     


    %Equations
    gradceq(4.*i+1,4.*i+1) = -1;
    gradceq(4.*i+1,4.*(i-1)+1) = 1;
    
    gradceq(4.*i+1,4.*(i-1)+3) = Inc_time.*dv;
    gradceq(4.*i+1,4.*(i-1)+4) = Inc_time.*dgamma;

    gradceq(4.*i+1,4.*h+1) = dt_f;

    %C2
    %Diferentials
    dv = sin(gamma_ant);
    dgamma = v_ant.*cos(gamma_ant);
    dt_f = (1./h).*v_ant.*sin(gamma_ant);

    %Equations
    gradceq(4.*i+2,4.*i+2) = -1;
    gradceq(4.*i+2,4.*(i-1)+2) = 1;
    
    gradceq(4.*i+2,4.*(i-1)+3) = Inc_time.*dv;
    gradceq(4.*i+2,4.*(i-1)+4) = Inc_time.*dgamma;

    gradceq(4.*i+2,4.*h+1) = dt_f;
    
    %C3
    
    dgamma = -g.*cos(gamma_ant);
    dt_f = -(1./h).*g.*sin(gamma_ant);


    gradceq(4.*i+3,4.*i+3) = -1;
    gradceq(4.*i+3,4.*(i-1)+3) = 1;
    
    
    gradceq(4.*i+3,4.*(i-1)+4) = Inc_time.*dgamma;

    gradceq(4.*i+3,4.*h+1) = dt_f;

    %C4
    dv = (g./(v_ant.^2)).*cos(gamma_ant);
    dgamma = (g./v_ant).*sin(gamma_ant);
    %dgamma_gamma = (1-(g./v_ant).*sin(gamma_ant))^(-1);
    dt_f = -(1./h).*(g./v_ant).*cos(gamma_ant);

    gradceq(4.*i+4,4.*i+4) = -1;
    gradceq(4.*i+4,4.*(i-1)+4) = 1 + Inc_time.*dgamma;
    
    
    gradceq(4.*i+4,4.*(i-1)+3) = Inc_time.*dv;

    gradceq(4.*i+4,4.*h+1) = (1./h).*dt_f;

    
end
%Variables
    x_1_ant = x_1_0;
    x_2_ant = x_2_0;
    v_ant = v_0;
    gamma_ant = x(4.*h+2);

    %Time variables:
    t_f = x(4.*h+1);

    Inc_time = ((t_f-t_0)/h);
    
    %C1
    %Diferentials
    %dgamma_gamma = (1-(g./v_ant).*sin(gamma_ant))^(-1);
    dgamma = -v_ant.*sin(gamma_ant);
    dt_f = (1./h).*v_ant.*cos(gamma_ant);
    %Equations:
    gradceq(1,1) = -1;
    gradceq(1,4.*h+1) = dt_f;
    gradceq(1,4.*h+2) = Inc_time.*dgamma;

    %C2
    %Diferentials
   
    dgamma = v_ant.*cos(gamma_ant);
    dt_f = (1./h).*v_ant.*sin(gamma_ant);
    %Equations
    gradceq(2,2) = -1;
    gradceq(2,4.*h+1) = dt_f;
    gradceq(2,4.*h+2) = Inc_time.*dgamma;

    %C3
    %Diferetials
    dgamma = -g.*cos(gamma_ant);
    dt_f = -(1./h).*g.*sin(gamma_ant);
    %Equations:

    gradceq(3,3) = -1;
    gradceq(3,4.*h+1) = dt_f;
    gradceq(3,4.*h+2) = Inc_time.*dgamma;

    %C4
    %Diferetials
    dgamma = (g./v_ant).*sin(gamma_ant);
    dt_f = -(1./h).*(g./v_ant).*cos(gamma_ant);
    %Equations
    gradceq(4,4) = -1;
    gradceq(4,4.*h+1) = dt_f;      
    gradceq(4,4.*h+2) = Inc_time.*dgamma;

    %C_last
    gradceq(4.*h+1,4.*(h-1)+2) = -1;

%transpose temporary

gradceq = transpose(gradceq);
end
