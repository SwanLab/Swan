function Avoigt = convert2Voigt(A,type)

    dim = size(A,1);
    if dim == 2
        pairs = [1 1; 2 2; 1 2];
    elseif dim == 3
        pairs = [1 1; 2 2; 3 3; 2 3; 1 3; 1 2];
    else
        error('Tensor must be 2D or 3D.');
    end
    dimVoigt = size(pairs,1);
    
    switch type
        case 'Strain'
            Avoigt = zeros(dimVoigt,1);
            for m = 1:dimVoigt
                    i = pairs(m,1); j = pairs(m,2);
                    Avoigt(m) = A(i,j);
            end
            nonDiag = (dim+1):length(Avoigt);
            Avoigt(nonDiag) = 2*Avoigt(nonDiag);
            
        case 'Stress'
            Avoigt = zeros(dimVoigt,1);
            for m = 1:dimVoigt
                    i = pairs(m,1); j = pairs(m,2);
                    Avoigt(m) = A(i,j);
            end

        case 'Constitutive'
            Avoigt = zeros(dimVoigt);
            for m = 1:dimVoigt
                for n = 1:dimVoigt
                    i = pairs(m,1); j = pairs(m,2);
                    k = pairs(n,1); l = pairs(n,2);
                    Avoigt(m,n) = A(i,j,k,l);
                end
            end

        case 'Compliance'
            Avoigt = zeros(dimVoigt);
            for m = 1:dimVoigt
                for n = 1:dimVoigt
                    i = pairs(m,1); j = pairs(m,2);
                    k = pairs(n,1); l = pairs(n,2);
                    Avoigt(m,n) = A(i,j,k,l);
                end
            end
            Avoigt(dim+1:end,:) = 2*Avoigt(dim+1:end,:);
            Avoigt(:,dim+1:end) = 2*Avoigt(:,dim+1:end);
    end

end