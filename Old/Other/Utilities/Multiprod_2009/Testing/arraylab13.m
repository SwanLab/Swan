function c = arraylab13(a,b,d1,d2)
% This is the engine used in MULTIPROD 1.3 for these cases:
% PxQ IN A - Rx1 IN B
% PxQ IN A - RxS IN B (slowest)

ndimsA = ndims(a); % NOTE - Since trailing singletons are removed,
ndimsB = ndims(b); %        not always NDIMSB = NDIMSA
NsA = d2 - ndimsA; % Number of added trailing singletons
NsB = d2 - ndimsB;
sizA = [size(a) ones(1,NsA)];
sizB = [size(b) ones(1,NsB)];

% Performing products
if sizB(d2) == 1 %  PxQ IN A  -  Rx1 IN B
    % A * B
    c = mbyv(a, b, d1);
else %          PxQ IN A  -  RxS IN B (least efficient)
    p = sizA(d1);
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
        Cindices = Bindices;
        Cindices{d1} = 1:p;            
    % Building C
        for Ncol = 1:s
           Bindices{d2} = Ncol; Cindices{d2} = Ncol;
           c(Cindices{:}) = mbyv(a, b(Bindices{:}), d1);
        end
end


function c = mbyv(a, b, dim)
% NOTE: This function is part of MULTIPROD 1.3

% 1 - Transposing: Qx1 matrices in B become 1xQ matrices
order = [1:dim-1, dim+1, dim, dim+2:ndims(b)];
b = permute(b, order);

% 2 - Cloning B P times along its singleton dimension DIM.
%        Optimized code for B2 = REPMAT(B, [ONES(1,DIM-1), P]).
%        INDICES is a cell array containing vectorized indices.
P = size(a, dim);
siz = size(b);
siz = [siz ones(1,dim-length(siz))]; % Ones are added if DIM > NDIMS(B)
Nd = length(siz);
indices = cell(1,Nd); % preallocating
for d = 1 : Nd
    indices{d} = 1:siz(d);
end
indices{dim} = ones(1, P); % "Cloned" index for dimension DIM
b2 = b(indices{:});        % B2 has same size as A

% 3 - Performing dot products along dimension DIM+1
c = sum(a .* b2, dim+1);