function y = launch(gamma0,tf,num_cons,phisical_cons)
g = phisical_cons(1,1);
v_0 =phisical_cons(1,2);
x_1_0 =phisical_cons(1,3);
x_2_0 =phisical_cons(1,4);
t_0 =phisical_cons(1,5);
h = num_cons;

%First set:

    dx_1 = v_0.*cos(gamma0);
    dx_2 = v_0.*sin(gamma0);
    d_v = -g.*sin(gamma0);
    d_gamma = -(g./v_0).*cos(gamma0);

%Inc of time

    Inc_time = ((tf-t_0)./h);
    
    
x_1_new = x_1_0 + Inc_time.*dx_1;
x_2_new = x_2_0 + Inc_time.*dx_2;
v_new = v_0 + Inc_time.*d_v;
gamma_new = gamma0 + Inc_time.*d_gamma;

y = zeros(4,h);
y(1,1) = x_1_0;
y(2,1) = x_2_0;
y(3,1) = v_0;
y(4,1) = gamma0;

    for i = 1:(h-1)
        x_1_ant = x_1_new;
        x_2_ant = x_2_new;
        v_ant = v_new;
        gamma_ant = gamma_new;

        dx_1 = v_ant.*cos(gamma_ant);
        dx_2 = v_ant.*sin(gamma_ant);
        d_v = -g.*sin(gamma_ant);
        d_gamma = -(g./v_ant).*cos(gamma_ant);

        x_1_new = x_1_ant + Inc_time.*dx_1;
        x_2_new = x_2_ant + Inc_time.*dx_2;
        v_new = v_ant + Inc_time.*d_v;
        gamma_new = gamma_ant + Inc_time.*d_gamma;

        y(1,i+1) = x_1_new;
        y(2,i+1) = x_2_new;
        y(3,i+1) = v_new;
        y(4,i+1) = gamma_new;

    end
end

