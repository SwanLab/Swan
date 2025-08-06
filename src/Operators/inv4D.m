function Ainv = inv4D(A)
    N    = size(A,1);
    Amat = reshape(A, [N^2, N^2]);
    Ainv = inv(Amat);
    Ainv = reshape(Ainv, [N, N, N, N]);
end