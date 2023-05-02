A = load("A_1000.mat").A_1000;
b = load("b_1000.mat").b_1000;

s = 100;
k = 1000;
tol = 10^(-4);

n   = size(A,1);
k   = 5000;
s   = 100;
tol = 10^(-4);

e1    = zeros(n,1);
e1(1) = 1;
x0 = obj.xPrevIt;
r0 = b - A*x0;

if exist(obj.Yk) == 0

    beta1  = norm(r0);
    V(:,1) = r0/beta1;
    c0     = beta1*e1;

    [xNew, rNew, V, H] = regularMINRES.solve(A, b, k, tol);
    x_k = xNew;
    r_k = rNew;
    Pk = createPk(k,s);
    obj.Yk = V*Pk;
else
    [Q, R] = qr(A*obj.Yk);
    Ck = Q;
    Uk = obj.Yk/R;
    x_start = obj.xPrevIt + Uk*Ck'*r0;
    r_start = r0-Ck*Ck'*r0;

    V(:,1) = r_start/norm(r_start);
    newA = (I-Ck*Ck')*A;

    [xNew, rNew, V, H] = regularMINRES.solve(newA, b, k ,tol);
    B = Ck'*A*V;
    norms = sqrt(sum(U.^2,1));
    Dk = diag(1./norms);
    W = [Ck V];
    G = [Dk B; 0 H];



    Pk = createPk(k,s);
    obj.Yk = V*Pk;

    %compute D

    norms = sqrt(sum(U.^2,1));
    Dk = diag(1./norms);












end