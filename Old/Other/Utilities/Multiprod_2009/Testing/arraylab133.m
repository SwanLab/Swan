function c = arraylab133(a,b,d1,d2)
% Several adjustments to ARRAYLAB13:
%    1) Adjustment used in ARRAYLAB131 was not used here.
%    2) Nested statement used in ARRAYLAB132 was used here.
%    3) PERMUTE in subfunction MBYV was substituted with RESHAPE
%       (faster by one order of magnitude!).

    ndimsA = ndims(a); % NOTE - Since trailing singletons are removed,
    ndimsB = ndims(b); %        not always NDIMSB = NDIMSA
    NsA = d2 - ndimsA; % Number of added trailing singletons
    NsB = d2 - ndimsB;
    sizA = [size(a) ones(1,NsA)];
    sizB = [size(b) ones(1,NsB)];
    p = sizA(d1);
    r = sizB(d1);
    s = sizB(d2); 
    % Initializing C
        sizC = sizA; 
        sizC(d2) = s;
        c = zeros(sizC);
    % Vectorized indices for B and C
        Nd = length(sizB);
        Bindices = cell(1,Nd); % preallocating (cell array)
        for d = 1 : Nd
           Bindices{d} = 1:sizB(d);
        end
        B2size = sizB; B2size([d1 d2]) = [1 r];
        B2indices = Bindices; 
        % B2 will be cloned P times along its singleton dimension D1 (see MBYV).
        B2indices([d1 d2]) = [{ones(1, p)} Bindices(d1)];  % "Cloned" index

    if sizB(d2) == 1 %  PxQ IN A  -  Rx1 IN B
        % A * B
        c = mbyv(a, b, B2indices,B2size,d1,d2,p);

    else %              PxQ IN A  -  RxS IN B
        Cindices = Bindices;
        Cindices{d1} = 1:p;            
        % Building C
        for Ncol = 1:s
           Bindices{d2} = Ncol; Cindices{d2} = Ncol;
           c(Cindices{:}) = mbyv(a, b(Bindices{:}), B2indices,B2size,d1,d2,p);
        end
    end

function c = mbyv(a, b2, indices, newsize, d1, d2, p)
% This is an adjustment to a subfunction used within MULTIPROD 1.3

% 1 - Transposing: Qx1 matrices in B become 1xQ matrices
b2 = reshape(b2, newsize);

% 3 - Performing dot products along dimension DIM+1
%      % NOTE: b(indices{:}) has same size as A
%      % NOTE: This nested statement is much faster than two separate ones.
c = sum(a .* b2(indices{:}), d2);