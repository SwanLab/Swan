function mat = IxI(n)
    I   = eye(n);
    mat = tensorprod(I,I);
end