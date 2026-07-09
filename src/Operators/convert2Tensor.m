function A = convert2Tensor(Avoigt,type)
    
    dimVoigt = size(Avoigt,1);
    if dimVoigt == 3
        dim = 2;
        pairs = [1 3; 3 2];
    elseif dimVoigt == 6
        dim = 3;
        pairs = [1 6 5; 6 2 4; 5 4 3];
    else
        error('Tensor must be 2D or 3D.');
    end

    
    switch type
        case 'Strain'
            A = zeros(dim,dim);
            for m=1:dim
                for n=1:dim
                    i = pairs(m,n);
                    A(m,n) = Avoigt(i);
                end
            end
            voigtCorrectionMatrix = 0.5*eye(size(A)) + 0.5;
            A = A.*voigtCorrectionMatrix;

        case 'Stress'
            A = zeros(dim,dim);
            for m=1:dim
                for n=1:dim
                    i = pairs(m,n);
                    A(m,n) = Avoigt(i);
                end
            end

        case 'Constitutive'
            A = zeros(dim,dim,dim,dim);
            for m=1:dim
                for n=1:dim
                    for o=1:dim
                        for p=1:dim
                            i = pairs(o,p); j=pairs(m,n);
                            A(m,n,o,p) = Avoigt(i,j);
                        end
                    end
                end
            end

        case 'Compliance'
            voigtCorrectionMatrix = ones(size(Avoigt));
            voigtCorrectionMatrix(dim+1:end,:) = 0.5*voigtCorrectionMatrix(dim+1:end,:);
            voigtCorrectionMatrix(:,dim+1:end) = 0.5*voigtCorrectionMatrix(:,dim+1:end);
            Avoigt = Avoigt.*voigtCorrectionMatrix;
            A = zeros(dim,dim,dim,dim);
            for m=1:dim
                for n=1:dim
                    for o=1:dim
                        for p=1:dim
                            i = pairs(o,p); j=pairs(m,n);
                            A(m,n,o,p) = Avoigt(i,j);
                        end
                    end
                end
            end
    end

end