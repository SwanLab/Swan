function [c, ceq,gradc,gradceq] = nlcon(x)

%Constants:
%g = phisical_cons(1,1);
%v_0 =phisical_cons(1,2);
%x_1_0 =phisical_cons(1,3);
%x_2_0 =phisical_cons(1,4);
%t_0 =phisical_cons(1,5);
%h = num_cons(1,1);
g = 9.81;
v_0 = 10;
x_1_0 = 0;
x_2_0 = 0;
t_0 = 0;
h = 300;

%inequality constrains
c = [];
% Initialize ceq as an array
    ceq = zeros(4.*h+1, 1);
%Constrains:
for i= 1:(h-1)
    ceq(4.*i+1,1) = x(4.*(i-1)+1) + ((x(4.*h+1)-t_0)/h).*x(4.*(i-1)+3).*cos(x(4.*(i-1)+4))-x(4.*i+1);
    ceq(4.*i+2,1) = x(4.*(i-1)+2) + ((x(4.*h+1)-t_0)/h).*x(4.*(i-1)+3).*sin(x(4.*(i-1)+4))-x(4.*i+2);
    ceq(4.*i+3,1) = x(4.*(i-1)+3) -((x(4.*h+1)-t_0)/h).*g.*sin(x(4.*(i-1)+4))-x(4.*i+3);
    ceq(4.*i+4,1) = x(4.*(i-1)+4) -((x(4.*h+1)-t_0)/h).*(g./x(4.*(i-1)+3)).*cos(x(4.*(i-1)+4))-x(4.*i+4);

end

%First set
ceq(1,1) = x_1_0 + ((x(4.*h+1)-t_0)/h).*v_0.*cos(x(4.*h+2))-x(1);
ceq(2,1) = x_2_0 + ((x(4.*h+1)-t_0)/h).*v_0.*sin(x(4.*h+2))-x(2);
ceq(3,1) = v_0 -((x(4.*h+1)-t_0)/h).*g.*sin(x(4.*h+2))-x(3);
ceq(4,1) = x(4.*h+2) -((x(4.*h+1)-t_0)/h).*(g./v_0).*cos(x(4.*h+2))-x(4);
%Second set

ceq(4.*h+1,1) = x(4.*(h-1)+2);

%Gradient calculation:
gradc = [];
gradceq = zeros(4.*h+1,4.*h+2); % For future code constrains in colums variables in rows
%Secction 1:
for i = 1:(h-1)
    %C1
    gradceq(4.*i+1,4.*i+1) = -1;
    gradceq(4.*i+1,4.*(i-1)+1) = 1;
    
    gradceq(4.*i+1,4.*(i-1)+3) = ((x(4.*h+1)-t_0)./h).*cos(4.*(i-1)+4);
    gradceq(4.*i+1,4.*(i-1)+4) = -((x(4.*h+1)-t_0)./h).*x(4.*(i-1)+3).*sin(4.*(i-1)+4);

    gradceq(4.*i+1,4.*h+1) = (1./h).*x(4.*(i-1)+3).*cos(x(4.*(i-1)+4));

    %C2
    gradceq(4.*i+2,4.*i+2) = -1;
    gradceq(4.*i+2,4.*(i-1)+2) = 1;
    
    gradceq(4.*i+2,4.*(i-1)+3) = ((x(4.*h+1)-t_0)./h).*sin(4.*(i-1)+4);
    gradceq(4.*i+2,4.*(i-1)+4) = ((x(4.*h+1)-t_0)./h).*x(4.*(i-1)+3).*cos(4.*(i-1)+4);

    gradceq(4.*i+2,4.*h+1) = (1./h).*x(4.*(i-1)+3).*sin(x(4.*(i-1)+4));
    
    %C3
    gradceq(4.*i+3,4.*i+3) = -1;
    gradceq(4.*i+3,4.*(i-1)+3) = 1;
    
    
    gradceq(4.*i+3,4.*(i-1)+4) = -((x(4.*h+1)-t_0)./h).*g.*cos(4.*(i-1)+4);

    gradceq(4.*i+3,4.*h+1) = -(1./h).*g.*sin(x(4.*(i-1)+4));

    %C4

    gradceq(4.*i+4,4.*i+4) = -1;
    gradceq(4.*i+4,4.*(i-1)+4) = 1 + ((x(4.*h+1)-t_0)./h).*(g./x(4.*(i-1)+3)).*sin(x(4.*(i-1)+4));
    
    
    gradceq(4.*i+4,4.*(i-1)+3) = ((x(4.*h+1)-t_0)./h).*(g./((x(4.*(i-1)+3)).^2)).*cos(4.*(i-1)+4);

    gradceq(4.*i+4,4.*h+1) = -(1./h).*(g./x(4.*(i-1)+3)).*cos(x(4.*(i-1)+4));

    
end
% C1
gradceq(1,1) = -1;
gradceq(1,4.*h+1) = (1./h).*v_0.*cos(x(4.*h+2));
gradceq(1,4.*h+2) = -((x(4.*h+1)-t_0)./h).*v_0.*sin(x(4.*h+2));

%C2
gradceq(2,2) = -1;
gradceq(2,4.*h+1) = (1./h).*v_0.*sin(x(4.*h+2));
gradceq(2,4.*h+2) = ((x(4.*h+1)-t_0)./h).*v_0.*cos(x(4.*h+2));

%C3
gradceq(3,3) = -1;
gradceq(3,4.*h+1) = -(1./h).*g.*sin(x(4.*h+2));
gradceq(3,4.*h+2) = -((x(4.*h+1)-t_0)./h).*g.*cos(x(4.*h+2));

%C4
gradceq(4,4) = -1;
gradceq(4,4.*h+1) = -(1./h).*(g./v_0).*cos(x(4.*h+2));      
gradceq(4,4.*h+2) = ((x(4.*h+1)-t_0)./h).*(g./v_0).*sin(x(4.*h+2));

%C_last
gradceq(4.*h+1,4.*(h-1)+2) = -1;

%transpose temporary

gradceq = transpose(gradceq);
end
