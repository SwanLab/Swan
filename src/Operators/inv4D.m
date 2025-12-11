function Ainv = inv4D(A)
    % INV4D   Invert a 4th-order 2D or 3D tensor using Voigt mapping
    % Input:  A(dim,dim,dim,dim)
    % Output: Ainv(dim,dim,dim,dim) - symmetrized inverse
    
    dim = size(A,1);
    
    % ----- Define Voigt mapping -----
    if dim == 2
        pairs = [1 1; 2 2; 1 2];
    elseif dim == 3
        pairs = [1 1; 2 2; 3 3; 2 3; 1 3; 1 2];
    else
        error('Tensor must be 2D or 3D.');
    end
    dimVoigt = size(pairs,1);
    
    % ----- Convert 4th-order tensor to Voigt matrix -----
    Avoigt = zeros(dimVoigt);
    for m = 1:dimVoigt
        for n = 1:dimVoigt
            i = pairs(m,1); j = pairs(m,2);
            k = pairs(n,1); l = pairs(n,2);
            Avoigt(m,n) = A(i,j,k,l);
        end
    end
    
    % ----- Invert Voigt matrix -----
    Avoigt_inv = inv(Avoigt);
    
    % ----- Convert back to 4th-order tensor -----
    Ainv = zeros(dim,dim,dim,dim);
    for m = 1:dimVoigt
        for n = 1:dimVoigt
            i = pairs(m,1); j = pairs(m,2);
            k = pairs(n,1); l = pairs(n,2);
    
            % Fill all minor and major symmetry positions
            vals = Avoigt_inv(m,n);
            Ainv(i,j,k,l) = vals;
            Ainv(j,i,k,l) = vals;
            Ainv(i,j,l,k) = vals;
            Ainv(j,i,l,k) = vals;
            Ainv(k,l,i,j) = vals;
            Ainv(l,k,i,j) = vals;
            Ainv(k,l,j,i) = vals;
            Ainv(l,k,j,i) = vals;
        end
    end
    
    for i=1:dim
        for j=1:dim
            if i ~= j
                Ainv(i,j,:,:) = Ainv(i,j,:,:)/4;
            end
        end
    end

end