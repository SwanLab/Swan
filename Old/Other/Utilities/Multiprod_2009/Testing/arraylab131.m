function c = arraylab131(a,b,d1,d2)
% Slight adjustment to ARRAYLAB13.
% I just used PERMUTE on the whole array, only once, rather than for each
% column of the RxS matrices. As a consequence, a call to subfunction MBYV
% is no longer necessary.

% Transposing: RxS matrices in B become SxR matrices
ndimsB = ndims(b);
order = [1:d1-1, d2, d1, d2+1:ndimsB];
b = permute(b, order);

% Checking
ndimsA = ndims(a); % NOTE - Since trailing singletons are removed,
ndimsB = ndims(b); %        not always NDIMSB = NDIMSA
NsA = d2 - ndimsA; % Number of added trailing singletons
NsB = d2 - ndimsB;
sizA = [size(a) ones(1,NsA)];
sizB = [size(b) ones(1,NsB)];

% Performing products
    p = sizA(d1); % A is PxQ
    s = sizB(d1); % B is SxR
    % Initializing C (PxS)
        sizC = sizA; 
        sizC(d2) = s;
        c = zeros(sizC);
    % Vectorized indices for B and C
        Nd = length(sizB);
        Bindices = cell(1,Nd); % preallocating (cell array)
        for d = 1 : Nd
           Bindices{d} = 1:sizB(d);
        end

    if sizB(d1) == 1 %  PxQ IN A  -  Rx1 IN B (after PERMUTing: 1xR in B)
        % A * B
        Bindices{d1} = ones(1, p); % "Cloned" index
        b2 = b(Bindices{:}); % B2 has same size as A
        c = sum(a .* b2, d2);
 
    else %              PxQ IN A  -  RxS IN B (after PERMUTing: SxR in B)

        Cindices = Bindices;
        Cindices{d1} = 1:p;            
        Cindices{d2} = Bindices{d1}; 
        % Building C
            for Ncol = 1:s
               % Cloning a part of B identified by index Ncol for dim. d1 
               Bindices{d1} = Ncol * ones(1, p); % "Cloned" index
               Cindices{d2} = Ncol;
               b2 = b(Bindices{:}); % B2 has same size as A

               % Performing dot products along dimension d2
               c(Cindices{:}) = sum(a .* b2, d2);
            end
    end

