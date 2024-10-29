% theta = sym('th');
% alpha = sym('alpha');
% A = [cosd(theta).^2 0; 0 ((sind(theta).^2)/alpha)];
% disp(inv(A));

theta = sym('th');
alpha = sym('alpha');
epsilon = sym('eps');

R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
a = [1/epsilon^2 0; 0 1/(epsilon * alpha)^2];

Result = R * a * R' ; 

disp(Result);
