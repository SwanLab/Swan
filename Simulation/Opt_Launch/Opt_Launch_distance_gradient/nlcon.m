function [c,ceq] = nlcon(x,num_cons,phisical_cons)
%Constants:
g = phisical_cons(1,1);
v_0 =phisical_cons(1,2);
x_1_0 =phisical_cons(1,3);
x_2_0 =phisical_cons(1,4);
t_0 =phisical_cons(1,5);
h = num_cons;

%variables
gamma0 = x(1);
tf = x(2);

y = launch(gamma0,tf,num_cons,phisical_cons);

%Inequality constrains:
c=[];
ceq = y(2,end)-0;
%gradc=[];
%gradceq=[];
y(2,end)

end

