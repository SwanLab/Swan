function [c, ceq,gradc,gradceq] = nlcon(x)

%Constants:
g = 9.81;
v_0 = 10;
x_1_0 = 0;
x_2_0 = 0;
t_0 = 0;
h = 3;

%inequality constrains
c = [];
% Initialize ceq as an array
    ceq = zeros(13, 1);
%Constrains:
%First set
ceq(1,1) = x_1_0 + ((x(13)-t_0)/h).*v_0.*cos(x(14))-x(1);
ceq(2,1) = x_2_0 + ((x(13)-t_0)/h).*v_0.*sin(x(14))-x(2);
ceq(3,1) = v_0 + -((x(13)-t_0)/h).*g.*sin(x(14))-x(3);
ceq(4,1) = x(14) + -((x(13)-t_0)/h).*(g./v_0).*cos(x(14))-x(4);
%Second set
ceq(5,1) = x(1) + ((x(13)-t_0)/h).*x(3).*cos(x(4))-x(5);
ceq(6,1) = x(2) + ((x(13)-t_0)/h).*x(3).*sin(x(4))-x(6);
ceq(7,1) = x(3) + -((x(13)-t_0)/h).*g.*sin(x(4))-x(7);
ceq(8,1) = x(4) + -((x(13)-t_0)/h).*(g./x(3)).*cos(x(4))-x(8);
%Third set
ceq(9,1) = x(5) + ((x(13)-t_0)/h).*x(7).*cos(x(8))-x(9);
ceq(10,1) = x(6) + ((x(13)-t_0)/h).*x(7).*sin(x(8))-x(10);
ceq(11,1) = x(7) + -((x(13)-t_0)/h).*g.*sin(x(8))-x(11);
ceq(12,1) = x(8) + -((x(13)-t_0)/h).*(g./x(7)).*cos(x(8))-x(12);
%Forth set
ceq(13,1) = x(10);

%Gradient calculation:
gradc = [];
gradceq = zeros(13,14); % For future code constrains in colums variables in rows
%Secction 1:

% C1
gradceq(1,1) = -1;
gradceq(1,13) = (1./h).*v_0.*cos(x(14));
gradceq(1,14) = -((x(13)-t_0)./h).*v_0.*sin(x(14));

%C2
gradceq(2,2) = -1;
gradceq(2,13) = (1./h).*v_0.*sin(x(14));
gradceq(2,14) = ((x(13)-t_0)./h).*v_0.*cos(x(14));

%C3
gradceq(3,3) = -1;
gradceq(3,13) = -(1./h).*g.*sin(x(14));
gradceq(3,14) = -((x(13)-t_0)./h).*g.*cos(x(14));

%C4
gradceq(4,4) = -1;
gradceq(4,13) = -(1./h).*(g./v_0).*cos(x(14));      
gradceq(4,14) = ((x(13)-t_0)./h).*(g./v_0).*sin(x(14));

%C5
gradceq(5,1) = 1;
gradceq(5,3) = ((x(13)-t_0)./h).*cos(x(4));
gradceq(5,4) = -((x(13)-t_0)./h).*x(3).*sin(x(4));

gradceq(5,5) = -1;
gradceq(5,13) = (1./h).*x(3).*cos(x(4));

%C6
gradceq(6,2) = 1;
gradceq(6,3) = ((x(13)-t_0)./h).*sin(x(4));
gradceq(6,4) = ((x(13)-t_0)./h).*x(3).*cos(x(4));

gradceq(6,6) = -1;
gradceq(6,13) = (1./h).*x(3).*sin(x(4));

%C7
gradceq(7,3) = 1;
gradceq(7,4) = -((x(13)-t_0)./h).*g.*cos(x(4));

gradceq(7,7) = -1;
gradceq(7,13) = -(1./h).*g.*sin(x(4));

%C8
gradceq(8,3) = ((x(13)-t_0)./h).*(g./x(3).^2).*cos(x(4));
gradceq(8,4) = 1+((x(13)-t_0)./h).*(g./x(3)).*sin(x(4));

gradceq(8,8) = -1;
gradceq(8,13) = -(1./h).*(g./x(3)).*sin(x(4));

%C9
gradceq(9,5) = 1;
gradceq(9,7) = ((x(13)-t_0)./h).*cos(x(8));
gradceq(9,8) = -((x(13)-t_0)./h).*x(7).*sin(x(8));

gradceq(9,9) = -1;
gradceq(9,13) = (1./h).*x(7).*cos(x(8));

%C10
gradceq(10,6) = 1;
gradceq(10,7) = ((x(13)-t_0)./h).*sin(x(8));
gradceq(10,8) = ((x(13)-t_0)./h).*x(7).*cos(x(8));

gradceq(10,10) = -1;
gradceq(10,13) = (1./h).*x(7).*sin(x(8));

%C11
gradceq(11,7) = 1;
gradceq(11,4) = -((x(13)-t_0)./h).*g.*cos(x(8));

gradceq(11,11) = -1;
gradceq(11,13) = -(1./h).*g.*sin(x(8));

%C12
gradceq(12,7) = -((x(13)-t_0)./h).*(g./x(7).^2).*cos(x(8));
gradceq(12,8) = 1+((x(13)-t_0)./h).*(g./x(7)).*sin(x(8));

gradceq(12,12) = -1;
gradceq(12,13) = -(1./h).*(g./x(7)).*sin(x(8));

%C13
gradceq(13,10) = -1;

%transpose temporary

gradceq = transpose(gradceq);
end
