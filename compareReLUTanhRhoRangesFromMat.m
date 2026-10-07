function compareReLUTanhRhoRangesFromMat()

clc;

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' ReLU vs tanh: dC/drho over restricted rho ranges\n');
fprintf('=============================================================\n');


%% ============================================================
% LOAD RELU NETWORK DIRECTLY FROM .MAT
% ============================================================

reluFile = 'HomogNN1.mat';

if ~exist(reluFile,'file')
    error('Could not find %s.',reluFile);
end

R = load(reluFile,'problem','data');

if ~isfield(R,'problem') || ~isfield(R,'data')
    error(['%s must contain variables problem and data.'],reluFile);
end

problemR = R.problem;
dataR    = R.data;


%% ============================================================
% LOAD TANH OBJECT
% ============================================================

tanhFile = 'BaselineCorrected_tanh_obj.mat';

if ~exist(tanhFile,'file')
    error('Could not find %s.',tanhFile);
end

T = load(tanhFile,'obj');
objT = T.obj;


%% ============================================================
% GRID / HOMOGENIZATION REFERENCE
% ============================================================

pB = objT.paramB(:)';
pR = objT.paramRho(:)';

C = objT.Chomog;

[BB,RR] = meshgrid(pB,pR);

fprintf('ReLU saved polynomial degree = %d\n',dataR.pol_deg);
fprintf('Grid = %d rho x %d b\n',numel(pR),numel(pB));


%% ============================================================
% CHECK SAVED RELU DOMAIN
% ============================================================

tol = 1e-10;

if abs(dataR.b_min-min(pB)) > tol || ...
   abs(dataR.b_max-max(pB)) > tol || ...
   abs(dataR.rho_min-min(pR)) > 1e-6 || ...
   abs(dataR.rho_max-max(pR)) > 1e-6

    warning(['Saved ReLU domain does not exactly match the tanh ', ...
             'homogenization grid. Check that HomogNN1.mat is ', ...
             'really the intended ReLU baseline.']);
end


%% ============================================================
% EVALUATE SAVED RELU dC/drho
%
% Returned order:
%
% 1 -> C1111
% 2 -> C2222
% 3 -> C1122
% 4 -> C1212
% 5 -> C1112
% 6 -> C2212
% ============================================================

dCdrhoR = evaluateSavedNetworkRhoDerivative( ...
    problemR,dataR,BB(:),RR(:));

nRho = numel(pR);
nB   = numel(pB);

dCdrhoR = reshape(dCdrhoR,nRho,nB,6);


%% ============================================================
% COMPONENTS
% ============================================================

ids = { ...
    [1 1 1 1], ...
    [1 1 2 2], ...
    [2 2 2 2], ...
    [1 2 1 2]};

names = { ...
    'C1111', ...
    'C1122', ...
    'C2222', ...
    'C1212'};

% Mapping to saved-network component order
compIndex = [1 3 2 4];

rhoCuts = [0.70 0.80 0.90 0.95 inf];


%% ============================================================
% COMPARISON
% ============================================================

for ic = 1:numel(ids)

    id = ids{ic};

    Cg = squeeze(C( ...
        id(1), ...
        id(2), ...
        id(3), ...
        id(4), ...
        :, ...
        :) );


    %% HOM reference derivative

    [~,dHOMdrho] = gradient(Cg,pB,pR);


    %% ReLU derivative reconstructed directly from HomogNN1.mat

    dRrho = dCdrhoR(:,:,compIndex(ic));


    %% tanh derivative from saved obj

    dT = objT.df{ ...
        id(1), ...
        id(2), ...
        id(3), ...
        id(4)}(BB,RR);

    dTrho = dT{2};


    fprintf('\n%s\n',names{ic});
    fprintf('-------------------------------------------------------------\n');
    fprintf('%10s %14s %14s %12s\n', ...
        'rho max','ReLU nRMSE','tanh nRMSE','tanh/ReLU');
    fprintf('-------------------------------------------------------------\n');


    for k = 1:numel(rhoCuts)

        rhoMax = rhoCuts(k);

        mask = false(size(BB));

        % Same interior region used in baselineDerivQuality
        mask(3:end-2,3:end-2) = true;

        if isfinite(rhoMax)
            mask(RR > rhoMax) = false;
        end


        ref = dHOMdrho(mask);
        rR  = dRrho(mask);
        rT  = dTrho(mask);


        nRelu = norm(rR-ref)/norm(ref);
        nTanh = norm(rT-ref)/norm(ref);


        if isfinite(rhoMax)
            rhoLabel = sprintf('<= %.2f',rhoMax);
        else
            rhoLabel = 'full';
        end


        fprintf('%10s %14.4e %14.4e %12.4f\n', ...
            rhoLabel, ...
            nRelu, ...
            nTanh, ...
            nTanh/nRelu);

    end
end


fprintf('\n');
fprintf('=============================================================\n');
fprintf(' Interpretation\n');
fprintf('=============================================================\n');
fprintf(' tanh/ReLU < 1  -> tanh has smaller derivative error\n');
fprintf(' tanh/ReLU > 1  -> ReLU has smaller derivative error\n');
fprintf('=============================================================\n');

end



% =============================================================
% SAVED NETWORK EVALUATOR
% =============================================================

function dCdrho = evaluateSavedNetworkRhoDerivative( ...
    problem,data,b,rho)


b = b(:);
rho = rho(:);

nPts = numel(b);


%% ============================================================
% NORMALIZE PHYSICAL VARIABLES
% ============================================================

bNorm = ...
    2*(b-data.b_min)/(data.b_max-data.b_min)-1;

rhoNorm = ...
    2*(rho-data.rho_min)/(data.rho_max-data.rho_min)-1;

Xbase = [bNorm,rhoNorm];


%% ============================================================
% POLYNOMIAL FEATURES
% ============================================================

[Xpoly,dXpolyDrhoNorm] = ...
    polynomialFeaturesAndRhoGradient( ...
    Xbase,data.pol_deg);


%% ============================================================
% STANDARDIZATION
% ============================================================

Xn = ...
    (Xpoly-data.muX)./data.stdX;


%% ============================================================
% PHYSICAL rho -> normalized rho SCALE
% ============================================================

scaleRho = ...
    2/(data.rho_max-data.rho_min);


dXdrho = ...
    dXpolyDrhoNorm*scaleRho;


% standardization derivative
dXdrho = ...
    dXdrho./data.stdX;


%% ============================================================
% NETWORK OUTPUT
% ============================================================

Yn = ...
    problem.computeOutputValues(Xn);

Y = ...
    Yn.*data.stdY + data.muY;


%% ============================================================
% NETWORK DIRECTIONAL DERIVATIVE wrt rho
% ============================================================

dX = ...
    reshape( ...
    dXdrho, ...
    [nPts,size(dXdrho,2),1]);


dYn = ...
    problem.computeDirectionalGradient(Xn,dX);


dYdrho = ...
    dYn(:,:,1).*data.stdY;


%% ============================================================
% CHOLESKY PARAMETERS
% ============================================================

L11 = exp(Y(:,1));
L21 = Y(:,2);

L22 = exp(Y(:,3));

L31 = Y(:,4);
L32 = Y(:,5);

L33 = exp(Y(:,6));


%% ============================================================
% THEIR DERIVATIVES
% ============================================================

dL11 = ...
    L11.*dYdrho(:,1);

dL21 = ...
    dYdrho(:,2);

dL22 = ...
    L22.*dYdrho(:,3);

dL31 = ...
    dYdrho(:,4);

dL32 = ...
    dYdrho(:,5);

dL33 = ...
    L33.*dYdrho(:,6);


%% ============================================================
% RECONSTRUCT dC/drho
%
% Same order as DamageHomogenizationFitter:
%
% 1 C1111
% 2 C2222
% 3 C1122
% 4 C1212
% 5 C1112
% 6 C2212
% ============================================================

dCdrho = zeros(nPts,6);


dCdrho(:,1) = ...
    2*L11.*dL11;


dCdrho(:,2) = ...
    2*L21.*dL21 ...
    + ...
    2*L22.*dL22;


dCdrho(:,3) = ...
    dL11.*L21 ...
    + ...
    L11.*dL21;


dCdrho(:,4) = ...
    2*L31.*dL31 ...
    + ...
    2*L32.*dL32 ...
    + ...
    2*L33.*dL33;


dCdrho(:,5) = ...
    dL11.*L31 ...
    + ...
    L11.*dL31;


dCdrho(:,6) = ...
    dL21.*L31 ...
    + ...
    L21.*dL31 ...
    + ...
    dL22.*L32 ...
    + ...
    L22.*dL32;

end



% =============================================================
% POLYNOMIAL FEATURES + DERIVATIVE wrt NORMALIZED rho
% =============================================================

function [Xpoly,dXdrho] = ...
    polynomialFeaturesAndRhoGradient(X,d)

N = size(X,1);

if d < 1
    Xpoly   = zeros(N,0);
    dXdrho  = zeros(N,0);
    return
end


Ecell = cell(d,1);

for g = 1:d
    Ecell{g} = generateExponents(2,g);
end

E = vertcat(Ecell{:});


e1 = E(:,1)';
e2 = E(:,2)';


x1 = X(:,1);
x2 = X(:,2);


Xpoly = ...
    x1.^e1 .* x2.^e2;


dXdrho = ...
    e2 ...
    .* ...
    x1.^e1 ...
    .* ...
    x2.^max(e2-1,0);

end



% =============================================================
% SAME EXPONENT ORDER AS DamageHomogenizationFitter
% =============================================================

function E = generateExponents(n,g)

if n == 1
    E = g;
    return
end

E = [];

for k = 0:g

    S = ...
        generateExponents(n-1,g-k);

    E = ...
        [E;
         k*ones(size(S,1),1),S];

end

end