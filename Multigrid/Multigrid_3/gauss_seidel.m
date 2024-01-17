function [u, res] = gauss_seidel(A, b, u, niter)

D = diag(A);
D = diag(D,0);
L = -(tril(A) - D);
U = -(triu(A) - D);
r_initial = norm(A*u - b); 
res = zeros(niter+1,1);
res(1,1) = 1;

for i = 1:niter+1
    u = (D - L) \ (U*u + b);   
    %r_curr = norm(A*u - b)/r_initial; 
    r_cur = norm(A*u - b,'inf');
    res(i+1) = r_cur; 
end

end
