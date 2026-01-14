dim = 2;
n = 1e5;
A = randn(dim,dim,n); %+1i*randn(3,3,n);

if dim==2
    tic
    D1 = eig2(A);
    t1=toc; % Elapsed time is 0.014055 seconds.
else
    tic
    D1 = eig3(A);
    t1=toc; % Elapsed time is 0.057625 seconds.
end

tic
D2 = zeros(dim,n);
for k=1:n
    D2(:,k) = eig(A(:,:,k));
end
t2=toc; % Elapsed time is 0.546626 (dim=2), 0.928497 seconds (dim=3).

fprintf('eig2/eig3 timing = %f [s] (n=%d)\n', t1, n);
fprintf('      eig timing = %f [s] (n=%d)\n', t2, n);

D1(:,1:10)
D2(:,1:10)