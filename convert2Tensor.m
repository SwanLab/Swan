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

        case 'Material'
            A = zeros(dim,dim,dim,dim);
            for m=1:dim
                for n=1:dim
                    for o=1:dim
                        for p=1:dim
                            i = pairs(m,n); j=pairs(o,p);
                            A(m,n,o,p) = Avoigt(i,j);
                        end
                    end
                end
            end

    end

    % A = zeros(dim,dim,dim,dim);
    % for m = 1:dimVoigt
    %     for n = 1:dimVoigt
    %         i = pairs(m,1); j = pairs(m,2);
    %         k = pairs(n,1); l = pairs(n,2);
    % 
    %         % Fill all minor and major symmetry positions
    %         vals = AvoigtInv(m,n);
    %         A(i,j,k,l) = vals;
    %         A(j,i,k,l) = vals;
    %         A(i,j,l,k) = vals;
    %         A(j,i,l,k) = vals;
    %         A(k,l,i,j) = vals;
    %         A(l,k,i,j) = vals;
    %         A(k,l,j,i) = vals;
    %         A(l,k,j,i) = vals;
    %     end
    % end
    % 
    % for i=1:dim
    %     for j=1:dim
    %         if i ~= j
    %             A(i,j,:,:) = A(i,j,:,:)/4;
    %         end
    %     end
    % end

end