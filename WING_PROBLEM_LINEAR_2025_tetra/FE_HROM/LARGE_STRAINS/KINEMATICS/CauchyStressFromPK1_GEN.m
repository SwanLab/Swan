function CauchyStress = CauchyStressFromPK1_GEN(PoneST,FgradST,detF,ndim)
%--------------------------------------------------------------------------
% Generalized version of CauchyStressFromPK1.m: accepts multiple columns in PoneST and FgradST
% Computes Cauchy stress tensor (Voigt notation) from PK1 stress and F
% Created by ChatGPT4, from  CauchyStressFromPK1.m. PROMPT: Here I have a function in which the inputs are assumed to be vectors.
% Generalize this function that it can handle matrices, and so the operations are performed across columns as well
%--------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end

if ndim == 2
    nF = 4;
    nstrain = 3;
elseif ndim == 3
    nF = 9;
    nstrain = 6;
else
    error('ndim must be 2 or 3');
end

% Handle multiple columns (each column = different scenario)
nCol = size(PoneST, 2);
nelem_ngaus = size(FgradST, 1) / nF;

% Create index maps
FROWS = cell(1, nF);
for i = 1:nF
    FROWS{i} = i:nF:nelem_ngaus*nF;
end

CROWS = cell(1, nstrain);
for i = 1:nstrain
    CROWS{i} = i:nstrain:nelem_ngaus*nstrain;
end

% Expand detF to match matrix shape if necessary
if isvector(detF)
    if size(detF,2) == 1 && nCol > 1
        detF = repmat(detF, 1, nCol);
    end
end

% Allocate result
CauchyStress = zeros(nelem_ngaus*nstrain, nCol);

% Apply expressions
if ndim == 2
    CauchyStress(CROWS{1},:) = (FgradST(FROWS{1},:) .* PoneST(FROWS{1},:) + ...
                                FgradST(FROWS{4},:) .* PoneST(FROWS{4},:)) ./ detF;
    
    CauchyStress(CROWS{2},:) = (FgradST(FROWS{2},:) .* PoneST(FROWS{2},:) + ...
                                FgradST(FROWS{3},:) .* PoneST(FROWS{3},:)) ./ detF;
    
    CauchyStress(CROWS{3},:) = (FgradST(FROWS{1},:) .* PoneST(FROWS{3},:) + ...
                                FgradST(FROWS{4},:) .* PoneST(FROWS{2},:)) ./ detF;

else % ndim == 3
    CauchyStress(CROWS{1},:) = (FgradST(FROWS{1},:) .* PoneST(FROWS{1},:) + ...
                                FgradST(FROWS{8},:) .* PoneST(FROWS{8},:) + ...
                                FgradST(FROWS{9},:) .* PoneST(FROWS{9},:)) ./ detF;

    CauchyStress(CROWS{2},:) = (FgradST(FROWS{2},:) .* PoneST(FROWS{2},:) + ...
                                FgradST(FROWS{6},:) .* PoneST(FROWS{6},:) + ...
                                FgradST(FROWS{7},:) .* PoneST(FROWS{7},:)) ./ detF;

    CauchyStress(CROWS{3},:) = (FgradST(FROWS{3},:) .* PoneST(FROWS{3},:) + ...
                                FgradST(FROWS{4},:) .* PoneST(FROWS{4},:) + ...
                                FgradST(FROWS{5},:) .* PoneST(FROWS{5},:)) ./ detF;

    CauchyStress(CROWS{4},:) = (FgradST(FROWS{2},:) .* PoneST(FROWS{4},:) + ...
                                FgradST(FROWS{7},:) .* PoneST(FROWS{3},:) + ...
                                FgradST(FROWS{6},:) .* PoneST(FROWS{5},:)) ./ detF;

    CauchyStress(CROWS{5},:) = (FgradST(FROWS{1},:) .* PoneST(FROWS{5},:) + ...
                                FgradST(FROWS{8},:) .* PoneST(FROWS{3},:) + ...
                                FgradST(FROWS{9},:) .* PoneST(FROWS{4},:)) ./ detF;

    CauchyStress(CROWS{6},:) = (FgradST(FROWS{1},:) .* PoneST(FROWS{6},:) + ...
                                FgradST(FROWS{9},:) .* PoneST(FROWS{2},:) + ...
                                FgradST(FROWS{8},:) .* PoneST(FROWS{7},:)) ./ detF;
end

end
