function tryingThirdOrder
x = sym('x');
%f = exp(x) - (x+2);
%f = sqrt(x) -sin(x);
%f = exp(-x)*sin(x);
f = sin(x)/x;
df = diff(f);
ddf = diff(df);

f = matlabFunction(f);
df = matlabFunction(df);
ddf = matlabFunction(ddf);

%f = @(x)  exp(x) - (x+2);
%df = @(x) exp(x) - 1;
%ddf = @(x) exp(x);
x0 = -100000;%0.8*pi;
xk = x0;
resNR  = 1;
iter = 1;
while resNR(iter) > 1e-14 && iter < 100
   
    b = df(xk);
    c = f(xk);
    incX_NR = -c/b;   
    xk = xk + incX_NR;
    iter = iter + 1;
    resNR(iter) = abs(f(xk));    
end

xk = x0;
resC = 1;
iter  = 1;
while resC(iter) > 1e-14 && iter < 100
    a = 0.5*ddf(xk);
    b = df(xk);
    c = f(xk);
    incX_C1  = (-b+sqrt(b*b-4*a*c))/(2*a);
    incX_C2  = (-b-sqrt(b*b-4*a*c))/(2*a);
    incX = min(incX_C1,incX_C2);
    xk = xk + incX;
    iter = iter + 1;    
    resC(iter) = abs(f(xk));
end

xk = x0;
resG = 1;
iter  = 1;
while resG(iter) > 1e-14   && iter < 100
    c = f(xk);
    incX_G = -c;   
    xk = xk + 0.1*incX_G;
    iter = iter + 1;
    resG(iter) = abs(f(xk)); 
end

figure()
semilogy(resC)
hold on
semilogy(resNR)
semilogy(resG)
legend('Cubic', 'NewtonRaphson','Gradient')


end
