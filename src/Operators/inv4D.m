function Ainv = inv4D(A)
    Amat = convertToVoigt(A);
    Ainv = inv(Amat);
    Ainv = convertFromVoigt(Ainv);
end

function Avoigt = convertToVoigt(A)
    dim = size(A,1);
    if dim == 2
        pairs = [1 1; 2 2; 1 2];
        dimVoigt = 3;
    elseif dim == 3
        pairs = [1 1; 2 2; 3 3; 2 3; 1 3; 1 2];
        dimVoigt = 6;
    else
        error("Matrix must be 2×2x2x2 (2D) or 3x3x3x3 (3D).");
    end
    Avoigt = zeros(dimVoigt);
    for m = 1:dimVoigt
        for n = 1:dimVoigt
            i=pairs(m,1); j=pairs(m,2);
            k=pairs(n,1); l=pairs(n,2);
            Avoigt(m,n) = A(i,j,k,l);
        end
    end
end

function A = convertFromVoigt(Avoigt)
    dimVoigt = size(Avoigt,1);
    if dimVoigt == 3
        dim = 2;
        pairs = [1 1; 2 2; 1 2];
    elseif dimVoigt == 6
        dim = 3;
        pairs = [1 1; 2 2; 3 3; 2 3; 1 3; 1 2];
    else
        error("Voigt matrix must be 3×3 (2D) or 6×6 (3D).");
    end

    A = zeros(dim,dim,dim,dim);
    for m = 1:dimVoigt
        for n = 1:dimVoigt
            i = pairs(m,1); j = pairs(m,2);
            k = pairs(n,1); l = pairs(n,2);
            A(i,j,k,l) = Avoigt(m,n);
        end
    end

    for i=1:dim
        for j=1:dim
            for k=1:dim
                for l=1:dim
                    A(i,j,k,l) = keepGreatestValue(A(k,l,i,j),A(i,j,k,l));
                    A(i,j,k,l) = keepGreatestValue(A(j,i,k,l),A(i,j,k,l));
                    A(i,j,k,l) = keepGreatestValue(A(i,j,l,k),A(i,j,k,l));
                end
            end
        end
    end

end

function c = keepGreatestValue(a,b)
   c = a * (abs(a) >= abs(b)) + b * (abs(b) > abs(a));
end