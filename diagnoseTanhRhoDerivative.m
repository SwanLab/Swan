function out = diagnoseTanhRhoDerivative()
% =============================================================
% diagnoseTanhRhoDerivative
%
% Spatial diagnostic of the rho derivative error for the
% corrected tanh surrogate.
%
% Compares:
%
%     dC_NN/drho
%
% against
%
%     dC_HOM/drho
%
% where the HOM derivative is obtained with MATLAB gradient()
% on the original homogenization grid.
%
% This script DOES NOT retrain anything.
%
% Output:
%   - global error metrics
%   - worst rho for each component
%   - worst b for each component
%   - location of maximum pointwise error
%   - percentage of total squared error concentrated in the
%     worst 5% of points
%   - heatmaps and line diagnostics
%   - saves results to TanhRhoDerivativeDiagnostic.mat
%
% =============================================================

clc;
close all;

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' TANH SURROGATE: SPATIAL DIAGNOSTIC OF dC/drho\n');
fprintf('=============================================================\n');


%% ============================================================
% LOAD SAVED TANH BASELINE
% ============================================================

fileName = ...
    'BaselineCorrected_tanh_obj.mat';

if ~exist(fileName,'file')

    error( ...
        'diagnoseTanhRhoDerivative:MissingFile', ...
        'Could not find %s.', ...
        fileName);

end

S = ...
    load(fileName,'obj');

obj = ...
    S.obj;


%% ============================================================
% DATA
% ============================================================

paramB = ...
    obj.paramB(:)';

paramRho = ...
    obj.paramRho(:)';

C = ...
    obj.Chomog;

dfun = ...
    obj.df;


nB = ...
    numel(paramB);

nRho = ...
    numel(paramRho);


[BB,RR] = ...
    meshgrid(paramB,paramRho);


fprintf('Grid: %d rho x %d b = %d points\n', ...
    nRho,nB,nRho*nB);

fprintf('b domain   : [%.8f, %.8f]\n', ...
    min(paramB),max(paramB));

fprintf('rho domain : [%.8f, %.8f]\n', ...
    min(paramRho),max(paramRho));


%% ============================================================
% SAME INTERIOR REGION USED BY baselineDerivQuality
%
% Exclude two grid layers at each boundary.
% ============================================================

interior = ...
    false(nRho,nB);

interior(3:end-2,3:end-2) = ...
    true;


rhoInterior = ...
    3:(nRho-2);

bInterior = ...
    3:(nB-2);


%% ============================================================
% COMPONENTS
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
% STORAGE
% ============================================================

out = ...
    struct();

out.fileName = ...
    fileName;

out.paramB = ...
    paramB;

out.paramRho = ...
    paramRho;

out.names = ...
    names;


out.results = ...
    struct([]);


%% ============================================================
% HEADER
% ============================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' GLOBAL / SPATIAL SUMMARY\n');
fprintf('=============================================================\n');

fprintf([ ...
    '%-7s %10s %9s %10s %10s %11s %11s\n'], ...
    'comp', ...
    'nRMSE', ...
    'sign', ...
    'rhoWorst', ...
    'bWorst', ...
    'top5energy', ...
    'max|err|');


%% ============================================================
% MAIN LOOP
% ============================================================

for ic = 1:nComp

    id = ...
        ids{ic};


    %% --------------------------------------------------------
    % HOMOGENIZATION COMPONENT
    % ---------------------------------------------------------

    Cg = ...
        squeeze( ...
        C( ...
        id(1), ...
        id(2), ...
        id(3), ...
        id(4), ...
        :, ...
        :) );


    %% --------------------------------------------------------
    % HOM DERIVATIVES
    %
    % Rows    -> rho
    % Columns -> b
    %
    % gradient returns:
    %
    %   first  output -> derivative along columns -> d/db
    %   second output -> derivative along rows    -> d/drho
    % ---------------------------------------------------------

    [~,dHOMdrho] = ...
        gradient(Cg,paramB,paramRho);


    %% --------------------------------------------------------
    % NN ANALYTICAL DERIVATIVE
    % ---------------------------------------------------------

    dNN = ...
        dfun{ ...
        id(1), ...
        id(2), ...
        id(3), ...
        id(4)}(BB,RR);


    if ~iscell(dNN) || numel(dNN) < 2

        error( ...
            'Unexpected dfun output for %s.', ...
            names{ic});

    end


    dNNdrho = ...
        dNN{2};


    %% --------------------------------------------------------
    % ERROR
    % ---------------------------------------------------------

    err = ...
        dNNdrho-dHOMdrho;


    refVec = ...
        dHOMdrho(interior);

    nnVec = ...
        dNNdrho(interior);

    errVec = ...
        err(interior);


    normRef = ...
        norm(refVec,2);


    if normRef == 0
        normRef = eps;
    end


    nRMSE = ...
        norm(errVec,2)/normRef;


    %% --------------------------------------------------------
    % SIGN MISMATCH
    %
    % Same criterion as baselineDerivQuality:
    %
    % |reference| > 5% max|reference|
    % ---------------------------------------------------------

    thresholdSign = ...
        0.05*max(abs(refVec));


    relevant = ...
        abs(refVec) > thresholdSign;


    if any(relevant)

        signMismatch = ...
            100*mean( ...
            sign(nnVec(relevant)) ...
            ~= ...
            sign(refVec(relevant)));

    else

        signMismatch = ...
            NaN;

    end


    %% --------------------------------------------------------
    % GLOBAL RMS REFERENCE SCALE
    %
    % Important:
    %
    % The per-rho and per-b errors below are normalized with ONE
    % global reference scale.
    %
    % We deliberately do NOT divide each row by its own local
    % derivative norm because a nearly-zero reference derivative
    % would artificially create huge ratios.
    % ---------------------------------------------------------

    globalRmsRef = ...
        sqrt(mean(refVec.^2));


    if globalRmsRef < eps
        globalRmsRef = eps;
    end


    %% --------------------------------------------------------
    % ERROR AS FUNCTION OF rho
    %
    % Each value measures RMS error across b at fixed rho,
    % normalized by the global HOM derivative RMS.
    % ---------------------------------------------------------

    errorVsRho = ...
        NaN(nRho,1);


    for iRho = rhoInterior

        eRow = ...
            err(iRho,bInterior);

        errorVsRho(iRho) = ...
            sqrt(mean(eRow.^2))/globalRmsRef;

    end


    [worstRhoScore,localIndex] = ...
        max(errorVsRho(rhoInterior));


    iWorstRho = ...
        rhoInterior(localIndex);


    rhoWorst = ...
        paramRho(iWorstRho);


    %% --------------------------------------------------------
    % ERROR AS FUNCTION OF b
    %
    % RMS across rho at fixed b.
    % ---------------------------------------------------------

    errorVsB = ...
        NaN(1,nB);


    for iB = bInterior

        eCol = ...
            err(rhoInterior,iB);

        errorVsB(iB) = ...
            sqrt(mean(eCol.^2))/globalRmsRef;

    end


    [worstBScore,localIndex] = ...
        max(errorVsB(bInterior));


    iWorstB = ...
        bInterior(localIndex);


    bWorst = ...
        paramB(iWorstB);


    %% --------------------------------------------------------
    % MAXIMUM POINTWISE ABSOLUTE ERROR
    % ---------------------------------------------------------

    absErr = ...
        abs(err);


    aux = ...
        absErr;

    aux(~interior) = ...
        -Inf;


    [maxAbsErr,linearIndex] = ...
        max(aux(:));


    [iMaxRho,iMaxB] = ...
        ind2sub(size(aux),linearIndex);


    bAtMax = ...
        paramB(iMaxB);


    rhoAtMax = ...
        paramRho(iMaxRho);


    refAtMax = ...
        dHOMdrho(iMaxRho,iMaxB);


    nnAtMax = ...
        dNNdrho(iMaxRho,iMaxB);


    %% --------------------------------------------------------
    % ERROR LOCALIZATION:
    %
    % What fraction of total squared error is contained in the
    % worst 5% of interior points?
    %
    % Interpretation:
    %
    % close to 100% -> error highly localized
    % small value   -> error distributed
    % ---------------------------------------------------------

    errEnergy = ...
        errVec.^2;


    errEnergy = ...
        sort(errEnergy,'descend');


    nTop = ...
        max(1,ceil(0.05*numel(errEnergy)));


    totalEnergy = ...
        sum(errEnergy);


    if totalEnergy > 0

        top5Energy = ...
            100*sum(errEnergy(1:nTop))/totalEnergy;

    else

        top5Energy = ...
            0;

    end


    %% --------------------------------------------------------
    % SOME MAGNITUDE INFORMATION
    % ---------------------------------------------------------

    rmsRef = ...
        sqrt(mean(refVec.^2));

    rmsNN = ...
        sqrt(mean(nnVec.^2));

    rmsErr = ...
        sqrt(mean(errVec.^2));


    maxRef = ...
        max(abs(refVec));

    maxNN = ...
        max(abs(nnVec));


    %% --------------------------------------------------------
    % SAVE COMPONENT RESULT
    % ---------------------------------------------------------

    r = ...
        struct();


    r.name = ...
        names{ic};

    r.C = ...
        Cg;

    r.dHOMdrho = ...
        dHOMdrho;

    r.dNNdrho = ...
        dNNdrho;

    r.error = ...
        err;

    r.absError = ...
        absErr;


    r.nRMSE = ...
        nRMSE;

    r.signMismatch = ...
        signMismatch;


    r.rmsReference = ...
        rmsRef;

    r.rmsNN = ...
        rmsNN;

    r.rmsError = ...
        rmsErr;


    r.maxReference = ...
        maxRef;

    r.maxNN = ...
        maxNN;


    r.errorVsRho = ...
        errorVsRho;

    r.errorVsB = ...
        errorVsB;


    r.rhoWorst = ...
        rhoWorst;

    r.rhoWorstIndex = ...
        iWorstRho;

    r.rhoWorstScore = ...
        worstRhoScore;


    r.bWorst = ...
        bWorst;

    r.bWorstIndex = ...
        iWorstB;

    r.bWorstScore = ...
        worstBScore;


    r.maxAbsError = ...
        maxAbsErr;

    r.bAtMaxError = ...
        bAtMax;

    r.rhoAtMaxError = ...
        rhoAtMax;

    r.referenceAtMaxError = ...
        refAtMax;

    r.nnAtMaxError = ...
        nnAtMax;


    r.top5EnergyPercent = ...
        top5Energy;


    if ic == 1
        out.results = r;
    else
        out.results(ic) = r;
    end


    %% --------------------------------------------------------
    % PRINT MAIN LINE
    % ---------------------------------------------------------

    fprintf( ...
        '%-7s %10.3e %8.2f%% %10.4f %10.4f %10.2f%% %11.3e\n', ...
        names{ic}, ...
        nRMSE, ...
        signMismatch, ...
        rhoWorst, ...
        bWorst, ...
        top5Energy, ...
        maxAbsErr);

end


%% ============================================================
% DETAILED WORST-POINT INFORMATION
% ============================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' WORST POINT OF EACH COMPONENT\n');
fprintf('=============================================================\n');


for ic = 1:nComp

    r = ...
        out.results(ic);


    fprintf('\n%s\n',r.name);

    fprintf( ...
        '  worst rho profile : rho = %.6f\n', ...
        r.rhoWorst);

    fprintf( ...
        '  worst b profile   : b   = %.6f\n', ...
        r.bWorst);

    fprintf( ...
        '  maximum |error|   : %.6e\n', ...
        r.maxAbsError);

    fprintf( ...
        '  location          : b = %.6f, rho = %.6f\n', ...
        r.bAtMaxError, ...
        r.rhoAtMaxError);

    fprintf( ...
        '  HOM derivative    : %.6e\n', ...
        r.referenceAtMaxError);

    fprintf( ...
        '  NN derivative     : %.6e\n', ...
        r.nnAtMaxError);

    fprintf( ...
        '  RMS HOM derivative: %.6e\n', ...
        r.rmsReference);

    fprintf( ...
        '  RMS NN derivative : %.6e\n', ...
        r.rmsNN);

    fprintf( ...
        '  RMS error         : %.6e\n', ...
        r.rmsError);

    fprintf( ...
        '  error energy in worst 5%% points: %.2f %%\n', ...
        r.top5EnergyPercent);

end


%% ============================================================
% COMPONENTS TO PLOT
%
% Main mechanical components.
% ============================================================

plotComp = ...
    [1 2 4 6];


%% ============================================================
% FIGURE 1
%
% HEATMAP:
%
% |dC_NN/drho - dC_HOM/drho|
%
% normalized by the GLOBAL RMS reference derivative.
%
% This is the most important visualization for locating errors.
% ============================================================

figure( ...
    'Name', ...
    'tanh dC-drho spatial error', ...
    'Color','w');


tiledlayout( ...
    2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


for ip = 1:numel(plotComp)

    ic = ...
        plotComp(ip);

    r = ...
        out.results(ic);


    normalizedAbsError = ...
        r.absError/max(r.rmsReference,eps);


    nexttile;


    imagesc( ...
        paramB, ...
        paramRho, ...
        normalizedAbsError);


    set(gca,'YDir','normal');

    hold on;


    plot( ...
        r.bAtMaxError, ...
        r.rhoAtMaxError, ...
        'wo', ...
        'MarkerSize',9, ...
        'LineWidth',1.5);


    xlabel('b');
    ylabel('\rho');

    title( ...
        sprintf( ...
        '%s: |error| / RMS(HOM)', ...
        r.name));


    colorbar;

    grid on;

end


%% ============================================================
% FIGURE 2
%
% ERROR RMS VS rho
%
% Tells us whether a specific density range is responsible.
% ============================================================

figure( ...
    'Name', ...
    'tanh dC-drho error versus rho', ...
    'Color','w');


tiledlayout( ...
    2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


for ip = 1:numel(plotComp)

    ic = ...
        plotComp(ip);

    r = ...
        out.results(ic);


    nexttile;


    plot( ...
        paramRho, ...
        r.errorVsRho, ...
        'LineWidth',1.5);


    hold on;


    plot( ...
        r.rhoWorst, ...
        r.rhoWorstScore, ...
        'o', ...
        'MarkerSize',8, ...
        'LineWidth',1.5);


    xlabel('\rho');

    ylabel( ...
        'RMS error / global RMS(HOM)');


    title(r.name);

    grid on;

end


%% ============================================================
% FIGURE 3
%
% ERROR RMS VS b
%
% Tells us whether a particular anisotropy region is responsible.
% ============================================================

figure( ...
    'Name', ...
    'tanh dC-drho error versus b', ...
    'Color','w');


tiledlayout( ...
    2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


for ip = 1:numel(plotComp)

    ic = ...
        plotComp(ip);

    r = ...
        out.results(ic);


    nexttile;


    plot( ...
        paramB, ...
        r.errorVsB, ...
        'LineWidth',1.5);


    hold on;


    plot( ...
        r.bWorst, ...
        r.bWorstScore, ...
        'o', ...
        'MarkerSize',8, ...
        'LineWidth',1.5);


    xlabel('b');

    ylabel( ...
        'RMS error / global RMS(HOM)');


    title(r.name);

    grid on;

end


%% ============================================================
% FIGURE 4
%
% At the worst b for each component, compare dC/drho as a
% function of rho.
%
% This is particularly useful for identifying spikes or
% excessively sharp transitions in rho.
% ============================================================

figure( ...
    'Name', ...
    'tanh dC-drho at worst b', ...
    'Color','w');


tiledlayout( ...
    2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


for ip = 1:numel(plotComp)

    ic = ...
        plotComp(ip);

    r = ...
        out.results(ic);


    iB = ...
        r.bWorstIndex;


    nexttile;


    plot( ...
        paramRho, ...
        r.dHOMdrho(:,iB), ...
        'LineWidth',1.6);


    hold on;


    plot( ...
        paramRho, ...
        r.dNNdrho(:,iB), ...
        '--', ...
        'LineWidth',1.6);


    xlabel('\rho');

    ylabel( ...
        '\partial C / \partial \rho');


    title( ...
        sprintf( ...
        '%s, b = %.3f', ...
        r.name, ...
        r.bWorst));


    legend( ...
        'HOM finite difference', ...
        'NN analytical', ...
        'Location','best');


    grid on;

end


%% ============================================================
% FIGURE 5
%
% At the worst rho for each component, compare dC/drho across b.
%
% This tells us whether the failure is global across b or
% concentrated at particular lattice distortions.
% ============================================================

figure( ...
    'Name', ...
    'tanh dC-drho at worst rho', ...
    'Color','w');


tiledlayout( ...
    2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


for ip = 1:numel(plotComp)

    ic = ...
        plotComp(ip);

    r = ...
        out.results(ic);


    iRho = ...
        r.rhoWorstIndex;


    nexttile;


    plot( ...
        paramB, ...
        r.dHOMdrho(iRho,:), ...
        'LineWidth',1.6);


    hold on;


    plot( ...
        paramB, ...
        r.dNNdrho(iRho,:), ...
        '--', ...
        'LineWidth',1.6);


    xlabel('b');

    ylabel( ...
        '\partial C / \partial \rho');


    title( ...
        sprintf( ...
        '%s, \\rho = %.3f', ...
        r.name, ...
        r.rhoWorst));


    legend( ...
        'HOM finite difference', ...
        'NN analytical', ...
        'Location','best');


    grid on;

end


%% ============================================================
% SAVE DIAGNOSTIC
% ============================================================

save( ...
    'TanhRhoDerivativeDiagnostic.mat', ...
    'out', ...
    '-v7.3');


fprintf('\n');
fprintf('=============================================================\n');
fprintf(' DIAGNOSTIC SAVED\n');
fprintf('=============================================================\n');

fprintf( ...
    'Saved as TanhRhoDerivativeDiagnostic.mat\n');


fprintf('\n');
fprintf('What to inspect first:\n');

fprintf('  1) GLOBAL / SPATIAL SUMMARY table\n');
fprintf('  2) rhoWorst values\n');
fprintf('  3) top5energy values\n');
fprintf('  4) heatmaps\n');
fprintf('  5) HOM vs NN curves at worst b and worst rho\n');

fprintf('=============================================================\n');

end