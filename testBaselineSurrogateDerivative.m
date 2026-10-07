function testBaselineSurrogateDerivative()

clc;
rng(1234);

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' TRAINED SURROGATE DERIVATIVE CONSISTENCY TEST\n');
fprintf('=============================================================\n');


%% ============================================================
% LOAD THE SAVED CORRECTED BASELINE
% ============================================================

fileName = ...
    'BaselineCorrected_ReLU_obj.mat';


if ~exist(fileName,'file')

    error( ...
        'testBaselineSurrogateDerivative:MissingFile', ...
        'Could not find %s.', ...
        fileName);

end


S = ...
    load(fileName,'obj');


obj = ...
    S.obj;


%% ============================================================
% EXTRACT THE TRAINED SURROGATE
% ============================================================

fun = ...
    obj.f;

dfun = ...
    obj.df;

paramB = ...
    obj.paramB(:);

paramRho = ...
    obj.paramRho(:);


bMin = ...
    min(paramB);

bMax = ...
    max(paramB);

rhoMin = ...
    min(paramRho);

rhoMax = ...
    max(paramRho);


rangeB = ...
    bMax-bMin;

rangeRho = ...
    rhoMax-rhoMin;


fprintf('b domain   : [%.8f, %.8f]\n', ...
    bMin,bMax);

fprintf('rho domain : [%.8f, %.8f]\n', ...
    rhoMin,rhoMax);


%% ============================================================
% TEST POINTS
%
% Use RANDOM OFF-GRID points.
%
% This is deliberate:
% we do not want to verify the derivative only at the original
% homogenization samples.
%
% We stay away from the domain boundaries so all central
% differences remain inside the training domain.
% ============================================================

nCheck = ...
    250;


marginFraction = ...
    0.05;


bLow = ...
    bMin + marginFraction*rangeB;

bHigh = ...
    bMax - marginFraction*rangeB;


rhoLow = ...
    rhoMin + marginFraction*rangeRho;

rhoHigh = ...
    rhoMax - marginFraction*rangeRho;


B = ...
    bLow ...
    + ...
    (bHigh-bLow)*rand(nCheck,1);


R = ...
    rhoLow ...
    + ...
    (rhoHigh-rhoLow)*rand(nCheck,1);


fprintf('\n');
fprintf('Random off-grid test points : %d\n',nCheck);
fprintf('Boundary margin             : %.1f %%\n', ...
    100*marginFraction);


%% ============================================================
% COMPONENTS
%
% Order is the same one we have been using in derivQuality.
% ============================================================

ids = { ...
    [1 1 1 1], ...
    [1 1 2 2], ...
    [1 1 1 2], ...
    [2 2 2 2], ...
    [2 2 1 2], ...
    [1 2 1 2]};


names = { ...
    'C1111', ...
    'C1122', ...
    'C1112', ...
    'C2222', ...
    'C2212', ...
    'C1212'};


nComp = ...
    numel(ids);


%% ============================================================
% FINITE-DIFFERENCE STEP SIZES
%
% h is relative to the physical range of each variable.
%
% Therefore:
%
%     h_b   = factor * (bMax-bMin)
%     h_rho = factor * (rhoMax-rhoMin)
%
% ============================================================

hFactors = ...
    [1e-2 ...
     1e-3 ...
     1e-4 ...
     1e-5 ...
     1e-6 ...
     1e-7 ...
     1e-8];


nH = ...
    numel(hFactors);


errB = ...
    zeros(nComp,nH);

errRho = ...
    zeros(nComp,nH);


absRmsB = ...
    zeros(nComp,nH);

absRmsRho = ...
    zeros(nComp,nH);


%% ============================================================
% ANALYTICAL DERIVATIVES
% ============================================================

analyticB = ...
    cell(nComp,1);

analyticRho = ...
    cell(nComp,1);


fprintf('\nEvaluating analytical derivatives...\n');


for ic = 1:nComp

    id = ...
        ids{ic};


    dfHandle = ...
        dfun{ ...
        id(1), ...
        id(2), ...
        id(3), ...
        id(4)};


    if isempty(dfHandle)

        error( ...
            'Empty dfun for component %s.', ...
            names{ic});

    end


    d = ...
        dfHandle(B,R);


    if ~iscell(d) || numel(d) < 2

        error( ...
            ['dfun for %s did not return the expected ', ...
             '{dC/db,dC/drho} cell.'], ...
            names{ic});

    end


    analyticB{ic} = ...
        d{1}(:);


    analyticRho{ic} = ...
        d{2}(:);


    if ...
            any(~isfinite(analyticB{ic})) ...
            || ...
            any(~isfinite(analyticRho{ic}))

        error( ...
            'Non-finite analytical derivative in %s.', ...
            names{ic});

    end

end


%% ============================================================
% FINITE-DIFFERENCE SWEEP
% ============================================================

fprintf('Evaluating finite differences...\n');


for ih = 1:nH

    factor = ...
        hFactors(ih);


    hB = ...
        factor*rangeB;


    hR = ...
        factor*rangeRho;


    for ic = 1:nComp

        id = ...
            ids{ic};


        fHandle = ...
            fun{ ...
            id(1), ...
            id(2), ...
            id(3), ...
            id(4)};


        %% ----------------------------------------------------
        % dC/db
        % -----------------------------------------------------

        CplusB = ...
            fHandle(B+hB,R);

        CminusB = ...
            fHandle(B-hB,R);


        fdB = ...
            (CplusB(:)-CminusB(:))/(2*hB);


        gB = ...
            analyticB{ic};


        diffB = ...
            gB-fdB;


        denomB = ...
            max( ...
                [norm(gB,2), ...
                 norm(fdB,2), ...
                 eps]);


        errB(ic,ih) = ...
            norm(diffB,2)/denomB;


        absRmsB(ic,ih) = ...
            sqrt(mean(diffB.^2));


        %% ----------------------------------------------------
        % dC/drho
        % -----------------------------------------------------

        CplusR = ...
            fHandle(B,R+hR);

        CminusR = ...
            fHandle(B,R-hR);


        fdR = ...
            (CplusR(:)-CminusR(:))/(2*hR);


        gR = ...
            analyticRho{ic};


        diffR = ...
            gR-fdR;


        denomR = ...
            max( ...
                [norm(gR,2), ...
                 norm(fdR,2), ...
                 eps]);


        errRho(ic,ih) = ...
            norm(diffR,2)/denomR;


        absRmsRho(ic,ih) = ...
            sqrt(mean(diffR.^2));

    end

end


%% ============================================================
% GLOBAL h SWEEP
%
% Report maximum and median error across the six tensor
% components for each h.
% ============================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' GLOBAL FINITE-DIFFERENCE SWEEP\n');
fprintf('=============================================================\n');

fprintf('\n');
fprintf('dC/db\n');

fprintf( ...
    '-----------------------------------------------------------------\n');

fprintf( ...
    ' factor         h_b          median rel.err      max rel.err\n');

fprintf( ...
    '-----------------------------------------------------------------\n');


for ih = 1:nH

    fprintf( ...
        '%8.1e   %12.4e      %12.4e      %12.4e\n', ...
        hFactors(ih), ...
        hFactors(ih)*rangeB, ...
        median(errB(:,ih)), ...
        max(errB(:,ih)));

end


fprintf( ...
    '-----------------------------------------------------------------\n');


fprintf('\n');
fprintf('dC/drho\n');

fprintf( ...
    '-----------------------------------------------------------------\n');

fprintf( ...
    ' factor         h_rho        median rel.err      max rel.err\n');

fprintf( ...
    '-----------------------------------------------------------------\n');


for ih = 1:nH

    fprintf( ...
        '%8.1e   %12.4e      %12.4e      %12.4e\n', ...
        hFactors(ih), ...
        hFactors(ih)*rangeRho, ...
        median(errRho(:,ih)), ...
        max(errRho(:,ih)));

end


fprintf( ...
    '-----------------------------------------------------------------\n');


%% ============================================================
% BEST h FOR EACH COMPONENT
% ============================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' BEST CONSISTENCY BY COMPONENT\n');
fprintf('=============================================================\n');


fprintf( ...
    ['%-7s  %12s  %10s  %12s   ', ...
     '%12s  %10s  %12s\n'], ...
    'comp', ...
    'rel.err db', ...
    'factor', ...
    'RMS abs db', ...
    'rel.err rho', ...
    'factor', ...
    'RMS abs rho');


bestErrB = ...
    zeros(nComp,1);

bestErrR = ...
    zeros(nComp,1);


bestFactorB = ...
    zeros(nComp,1);

bestFactorR = ...
    zeros(nComp,1);


for ic = 1:nComp

    [bestErrB(ic),iBestB] = ...
        min(errB(ic,:));


    [bestErrR(ic),iBestR] = ...
        min(errRho(ic,:));


    bestFactorB(ic) = ...
        hFactors(iBestB);


    bestFactorR(ic) = ...
        hFactors(iBestR);


    fprintf( ...
        ['%-7s  %12.4e  %10.1e  %12.4e   ', ...
         '%12.4e  %10.1e  %12.4e\n'], ...
        names{ic}, ...
        bestErrB(ic), ...
        bestFactorB(ic), ...
        absRmsB(ic,iBestB), ...
        bestErrR(ic), ...
        bestFactorR(ic), ...
        absRmsRho(ic,iBestR));

end


%% ============================================================
% WORST COMPONENT
% ============================================================

[maxBestB,iWorstB] = ...
    max(bestErrB);


[maxBestR,iWorstR] = ...
    max(bestErrR);


fprintf('\n');
fprintf('Worst best-case dC/db consistency:\n');

fprintf( ...
    '  %s : %.6e\n', ...
    names{iWorstB}, ...
    maxBestB);


fprintf('\n');
fprintf('Worst best-case dC/drho consistency:\n');

fprintf( ...
    '  %s : %.6e\n', ...
    names{iWorstR}, ...
    maxBestR);


%% ============================================================
% SIGN CHECK AT EACH COMPONENT'S BEST h
%
% Only compare signs where the FD derivative is significant:
%
%       |FD| > 5 %% max|FD|
%
% This avoids counting numerical sign flips around zero.
% ============================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' SIGN CONSISTENCY AT BEST h\n');
fprintf('=============================================================\n');


fprintf( ...
    '%-7s %12s %14s\n', ...
    'comp', ...
    'sign db', ...
    'sign drho');


for ic = 1:nComp

    id = ...
        ids{ic};


    fHandle = ...
        fun{ ...
        id(1), ...
        id(2), ...
        id(3), ...
        id(4)};


    %% b

    hB = ...
        bestFactorB(ic)*rangeB;


    fdB = ...
        ( ...
        fHandle(B+hB,R) ...
        - ...
        fHandle(B-hB,R) ...
        )/(2*hB);


    fdB = ...
        fdB(:);


    gB = ...
        analyticB{ic};


    thresholdB = ...
        0.05*max(abs(fdB));


    relevantB = ...
        abs(fdB) > thresholdB;


    if any(relevantB)

        signMismatchB = ...
            100*mean( ...
            sign(gB(relevantB)) ...
            ~= ...
            sign(fdB(relevantB)));

    else

        signMismatchB = ...
            NaN;

    end


    %% rho

    hR = ...
        bestFactorR(ic)*rangeRho;


    fdR = ...
        ( ...
        fHandle(B,R+hR) ...
        - ...
        fHandle(B,R-hR) ...
        )/(2*hR);


    fdR = ...
        fdR(:);


    gR = ...
        analyticRho{ic};


    thresholdR = ...
        0.05*max(abs(fdR));


    relevantR = ...
        abs(fdR) > thresholdR;


    if any(relevantR)

        signMismatchR = ...
            100*mean( ...
            sign(gR(relevantR)) ...
            ~= ...
            sign(fdR(relevantR)));

    else

        signMismatchR = ...
            NaN;

    end


    fprintf( ...
        '%-7s %10.2f %% %12.2f %%\n', ...
        names{ic}, ...
        signMismatchB, ...
        signMismatchR);

end


%% ============================================================
% FINAL DECISION
% ============================================================

overallError = ...
    max([bestErrB;bestErrR]);


fprintf('\n');
fprintf('=============================================================\n');
fprintf(' FINAL RESULT\n');
fprintf('=============================================================\n');

fprintf( ...
    'Maximum best-case relative error = %.6e\n', ...
    overallError);


if overallError < 1e-5

    fprintf('\nRESULT: PASS\n');

    fprintf( ...
        ['The analytical surrogate derivatives are ', ...
         'consistent with finite differences of the ', ...
         'trained surrogate.\n']);

elseif overallError < 1e-3

    fprintf('\nRESULT: ACCEPTABLE / INSPECT\n');

    fprintf( ...
        ['The derivative chain is broadly consistent, ', ...
         'but ReLU activation boundaries may affect ', ...
         'some sampled points.\n']);

else

    fprintf('\nRESULT: FAIL / INVESTIGATE\n');

    fprintf( ...
        ['The analytical derivative chain does not ', ...
         'match finite differences closely enough.\n']);

end


fprintf('=============================================================\n');

end