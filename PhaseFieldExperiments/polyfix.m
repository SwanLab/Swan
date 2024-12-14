function P = polyfix(xi,yi,x0,y0,m)

% P = polyfix(xi,yi,x0,y0,m)
% This function fits a polynomial passing thru a fixed point of (x0, y0)
% (xi, yi) are the data points to be fited.
% It will give P = [P1,P2,P3,...,Pm,Pm+1] corresponding to the polynomial of:
% y = P1 x^m + P2 x^(m-1) + ,..., + Pm x + Pm+1

a = zeros(m);
b = zeros(m,1);

%L = length(xi);
%N = 1:L;
for row = 1:m
    for col = 1:m
        a(row, col) = sum((xi-x0).^col.*(xi-x0).^row);
    end
    b(row) = sum((yi-y0).*(xi-x0).^row);
end
p  = (a\b)';

P = zeros(1, m+1);

for M = 0:m
    if M == 0
        K = 1:m;
    else
        K = M:m;
    end
    P(m+1-M) = sum(factorial(K)./factorial(K-M)./factorial(M).*p(K).*(-x0).^(K-M));
end
P(m+1) = P(m+1) + y0;

return