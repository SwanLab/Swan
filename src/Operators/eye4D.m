function mat = eye4D(n)
    dik = reshape(1:n, [n, 1, 1, 1]) == reshape(1:n, [1, 1, n, 1]);
    djl = reshape(1:n, [1, n, 1, 1]) == reshape(1:n, [1, 1, 1, n]);
    mat  = double(dik & djl);
end
 