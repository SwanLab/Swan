clear;

%coordinates of the triangle vertex (x1, y1, x2, y2, etc) (random)
coordElem = [ 2 2; 7 1; 5 4];

% Compute area of the element
Ae = 0.5*abs(((coordElem(2,1)*coordElem(3,2)-coordElem(3,1)*coordElem(2,2))+(coordElem(2,2)-coordElem(3,2))*coordElem(1,1)+(coordElem(3,1)-coordElem(2,1))*coordElem(1,2)));

% matrix of coefficients
coef = zeros(length(coordElem),3);

% coeficients of linear equations
for i = 1 : length(coordElem)

    if i == 1
        j = 2; k = 3;
    elseif i == 2
        j = 3; k = 1;
    elseif i == 3
        j = 1; k = 2;
    end

    coef(i,1) = (1/(2*Ae)) * ( coordElem(j,2) - coordElem(k,2) );
    coef(i,2) = (1/(2*Ae)) * ( coordElem(k,1) - coordElem(j,1) );
    coef(i,3) = (1/(2*Ae)) * ( coordElem(j,1) * coordElem(k,2) - coordElem(k,1) * coordElem(j,2) );

end

%AD using my code
x = ValDerForward(0,[1 0]);
y = ValDerForward(0,[0 1]);

N(1) = coef(1,1)*x + coef(1,2)*y + coef(1,3);
N(2) = coef(2,1)*x + coef(2,2)*y + coef(2,3);
N(3) = coef(3,1)*x + coef(3,2)*y + coef(3,3);

n = [N(1).double; N(2).double; N(3).double];

%Jacobian
J = transpose(n(:,2:end)) * coordElem;
