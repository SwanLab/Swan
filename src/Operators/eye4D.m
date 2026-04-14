function mat = eye4D(n)
    dik = reshape(1:n, [n, 1, 1, 1]) == reshape(1:n, [1, 1, n, 1]);
    djl = reshape(1:n, [1, n, 1, 1]) == reshape(1:n, [1, 1, 1, n]);
    dil = reshape(1:n, [n, 1, 1, 1]) == reshape(1:n, [1, 1, 1, n]);
    djk = reshape(1:n, [1, n, 1, 1]) == reshape(1:n, [1, 1, n, 1]);
    mat  = 0.5*(double(dik & djl) + double(dil & djk));
end
 