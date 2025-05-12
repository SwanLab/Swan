function mat = kronEye(n)
    I   = eye(n);
    mat = tensorprod(I,I);
end