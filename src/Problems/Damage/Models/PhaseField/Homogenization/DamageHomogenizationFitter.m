classdef DamageHomogenizationFitter < handle

    methods (Access = public, Static)

        function [fun, dfun, ddfun] = computePolynomial(degPoly, phi, C)
            obj = DamageHomogenizationFitter();
            fun = obj.computeFitting(degPoly, phi, C);
            [dfun, ddfun] = obj.computeDerivative(fun);
            [fun, dfun, ddfun] = obj.convertToHandle(fun, dfun, ddfun);
        end

        function [fun, dfun, ddfun] = computeNN(paramVectors, C, varargin)

            if ~iscell(paramVectors) || isempty(paramVectors)
                error('DamageHomogenizationFitter:InvalidParameters', ...
                    'paramVectors must be a non-empty cell array.');
            end

            s = struct();
            if nargin >= 3 && isstruct(varargin{1})
                s = varargin{1};
            end

            cfg = DamageHomogenizationFitter.createConfig(paramVectors, s);

            if cfg.retrain

                rng(cfg.seed, 'twister');

                data = DamageHomogenizationFitter.prepareData( ...
                    paramVectors, C, cfg);

                problem = DamageHomogenizationFitter.trainNetwork( ...
                    data, cfg);

                history = problem.getHistory();

                save(cfg.saveFile, ...
                    'problem', 'data', 'history', 'cfg', '-v7.3');

                save(cfg.historyFile, ...
                    'history', 'cfg', '-v7.3');

                fprintf('Network saved in %s\n', cfg.saveFile);

            else

                if ~exist(cfg.saveFile, 'file')
                    error('DamageHomogenizationFitter:MissingModel', ...
                        'Model file %s does not exist.', cfg.saveFile);
                end

                loaded = load(cfg.saveFile, ...
                    'problem', 'data', 'cfg');

                if ~isfield(loaded, 'problem') || ...
                   ~isfield(loaded, 'data') || ...
                   ~isfield(loaded, 'cfg')
                    error('DamageHomogenizationFitter:LegacyModel', ...
                        ['The saved model does not contain problem, data and cfg. ', ...
                         'Retrain it with the current fitter.']);
                end

                DamageHomogenizationFitter.validateLoadedModel( ...
                    paramVectors, cfg, loaded.data, loaded.cfg);

                problem = loaded.problem;
                data = loaded.data;

                fprintf('Network loaded from %s\n', cfg.saveFile);

            end

            [fun, dfun, ddfun] = ...
                DamageHomogenizationFitter.buildHandles(problem, data);

        end

    end

    methods (Access = private)

        function fun = computeFitting(~, degPoly, phi, C)

            phi = reshape(phi, length(phi), []);
            nStre = size(C, 1);
            fun = cell(2,2,2,2);

            for i = 1:nStre
                for j = 1:nStre
                    for k = 1:nStre
                        for l = 1:nStre

                            coeffs = polyfit( ...
                                phi, ...
                                squeeze(C(i,j,k,l,:)), ...
                                degPoly);

                            fun{i,j,k,l} = poly2sym(coeffs);

                            if isempty(symvar(fun{i,j,k,l}))
                                syms x
                                fun{i,j,k,l} = 1e-20*x.^9;
                            end

                        end
                    end
                end
            end

        end

        function [dfun, ddfun] = computeDerivative(~, fun)

            nStre = size(fun, 1);
            dfun = cell(2,2,2,2);
            ddfun = cell(2,2,2,2);

            for i = 1:nStre
                for j = 1:nStre
                    for k = 1:nStre
                        for l = 1:nStre
                            dfun{i,j,k,l} = diff(fun{i,j,k,l});
                            ddfun{i,j,k,l} = diff(dfun{i,j,k,l});
                        end
                    end
                end
            end

        end

        function [fun, dfun, ddfun] = convertToHandle(~, fun, dfun, ddfun)

            nStre = size(fun, 1);

            for i = 1:nStre
                for j = 1:nStre
                    for k = 1:nStre
                        for l = 1:nStre
                            fun{i,j,k,l} = matlabFunction(fun{i,j,k,l});
                            dfun{i,j,k,l} = matlabFunction(dfun{i,j,k,l});
                            ddfun{i,j,k,l} = matlabFunction(ddfun{i,j,k,l});
                        end
                    end
                end
            end

        end

    end

    methods (Access = private, Static)

        function cfg = createConfig(paramVectors, s)

            nVar = numel(paramVectors);

            cfg.retrain = false;
            cfg.seed = 0;

            cfg.parameterNames = arrayfun( ...
                @(k) sprintf('p%d', k), ...
                1:nVar, ...
                'UniformOutput', false);

            cfg.transforms = repmat({'identity'}, 1, nVar);

            cfg.featureMap = 'direct';
            cfg.pol_deg = 1;

            cfg.hiddenLayers = [50 100 200 100 50 30];
            cfg.HUtype = 'ReLU';
            cfg.OUtype = 'linear';

            cfg.optimizerType = 'SGD';
            cfg.maxEpochs = 500000;
            cfg.learningRate = 0.015;
            cfg.earlyStop = [];

            cfg.costType = 'L2';
            cfg.lambda = 0;

            cfg.saveFile = 'HomogNN.mat';
            cfg.historyFile = 'NNhistory.mat';
            cfg.referenceDataFile = '';

            cfg.allowExtrapolation = false;

            fields = fieldnames(s);

            for i = 1:numel(fields)
                name = fields{i};

                if ~isfield(cfg, name)
                    error('DamageHomogenizationFitter:UnknownOption', ...
                        'Unknown configuration option: %s.', name);
                end

                cfg.(name) = s.(name);
            end

            cfg.featureMap = lower(char(cfg.featureMap));
            cfg.optimizerType = char(cfg.optimizerType);
            cfg.HUtype = char(cfg.HUtype);
            cfg.OUtype = char(cfg.OUtype);
            cfg.costType = char(cfg.costType);
            cfg.saveFile = char(cfg.saveFile);
            cfg.historyFile = char(cfg.historyFile);
            cfg.referenceDataFile = char(cfg.referenceDataFile);

            if isstring(cfg.parameterNames)
                cfg.parameterNames = cellstr(cfg.parameterNames);
            end

            if isstring(cfg.transforms)
                cfg.transforms = cellstr(cfg.transforms);
            end

            if numel(cfg.parameterNames) ~= nVar
                error('DamageHomogenizationFitter:InvalidParameterNames', ...
                    'parameterNames must contain one name per variable.');
            end

            if numel(cfg.transforms) ~= nVar
                error('DamageHomogenizationFitter:InvalidTransforms', ...
                    'transforms must contain one transform per variable.');
            end

            for q = 1:nVar
                cfg.parameterNames{q} = char(cfg.parameterNames{q});
                cfg.transforms{q} = lower(char(cfg.transforms{q}));

                if ~ismember(cfg.transforms{q}, ...
                        {'identity', 'log', 'atanh'})
                    error('DamageHomogenizationFitter:InvalidTransform', ...
                        'Unsupported transform %s.', cfg.transforms{q});
                end
            end

            if ~ismember(cfg.featureMap, {'direct', 'polynomial'})
                error('DamageHomogenizationFitter:InvalidFeatureMap', ...
                    'featureMap must be direct or polynomial.');
            end

            if strcmp(cfg.featureMap, 'polynomial')
                if cfg.pol_deg < 1 || ...
                   abs(cfg.pol_deg - round(cfg.pol_deg)) > 0
                    error('DamageHomogenizationFitter:InvalidPolynomialDegree', ...
                        'pol_deg must be a positive integer.');
                end
            end

            if cfg.maxEpochs < 1
                error('DamageHomogenizationFitter:InvalidMaxEpochs', ...
                    'maxEpochs must be positive.');
            end

            if cfg.learningRate <= 0
                error('DamageHomogenizationFitter:InvalidLearningRate', ...
                    'learningRate must be positive.');
            end

            if cfg.lambda < 0
                error('DamageHomogenizationFitter:InvalidLambda', ...
                    'lambda must be non-negative.');
            end

        end

        function data = prepareData(paramVectors, C, cfg)

            nVar = numel(paramVectors);

            DamageHomogenizationFitter.validateParameterVectors( ...
                paramVectors);

            DamageHomogenizationFitter.validateTensorDimensions( ...
                paramVectors, C);

            grid = cell(1, nVar);
            axesColumn = cell(1, nVar);

            for q = 1:nVar
                axesColumn{q} = paramVectors{q}(:);
            end

            [grid{:}] = ndgrid(axesColumn{:});

            nPts = numel(grid{1});
            P = zeros(nPts, nVar);

            for q = 1:nVar
                P(:,q) = grid{q}(:);
            end

            Y = DamageHomogenizationFitter.computeCholeskyTargets(C, nVar);

            [U, ~] = DamageHomogenizationFitter.transformParameters( ...
                P, cfg.transforms);

            uMin = min(U, [], 1);
            uMax = max(U, [], 1);
            uRange = uMax - uMin;

            if any(uRange <= 0)
                error('DamageHomogenizationFitter:DegenerateDomain', ...
                    'Every parameter must vary over a non-zero interval.');
            end

            Z = 2*(U - uMin)./uRange - 1;

            [Xraw, ~] = ...
                DamageHomogenizationFitter.buildFeaturesAndDerivatives( ...
                Z, cfg.featureMap, cfg.pol_deg);

            [iTrain, iValidation] = ...
                DamageHomogenizationFitter.createSplit( ...
                P, cfg.referenceDataFile);

            muX = mean(Xraw(iTrain,:), 1);
            stdX = std(Xraw(iTrain,:), 0, 1);
            stdX(stdX == 0) = 1;

            muY = mean(Y(iTrain,:), 1);
            stdY = std(Y(iTrain,:), 0, 1);
            stdY(stdY == 0) = 1;

            Xn = (Xraw - muX)./stdX;
            Yn = (Y - muY)./stdY;

            data.Xtrain = Xn(iTrain,:);
            data.Ytrain = Yn(iTrain,:);

            data.Xvalidation = Xn(iValidation,:);
            data.Yvalidation = Yn(iValidation,:);

            data.Xtest = data.Xvalidation;
            data.Ytest = data.Yvalidation;

            data.iTrain = iTrain;
            data.iValidation = iValidation;

            data.parameterPoints = P;

            data.paramVectors = cell(1, nVar);
            for q = 1:nVar
                data.paramVectors{q} = paramVectors{q}(:)';
            end

            data.parameterNames = cfg.parameterNames;
            data.transforms = cfg.transforms;

            data.featureMap = cfg.featureMap;
            data.pol_deg = cfg.pol_deg;

            data.pMin = min(P, [], 1);
            data.pMax = max(P, [], 1);

            data.uMin = uMin;
            data.uMax = uMax;

            data.muX = muX;
            data.stdX = stdX;

            data.muY = muY;
            data.stdY = stdY;

            data.nVariables = nVar;
            data.nFeatures = size(Xn, 2);
            data.nLabels = 6;

            data.allowExtrapolation = cfg.allowExtrapolation;

            data.componentMap = { ...
                [1 1 1 1], ...
                [2 2 2 2], ...
                [1 1 2 2], ...
                [1 2 1 2], ...
                [1 1 1 2], ...
                [2 2 1 2]};

        end

        function validateParameterVectors(paramVectors)

            for q = 1:numel(paramVectors)

                p = paramVectors{q};

                if ~isnumeric(p) || isempty(p) || ~isvector(p)
                    error('DamageHomogenizationFitter:InvalidParameterVector', ...
                        'Each parameter vector must be a non-empty numeric vector.');
                end

                if any(~isfinite(p(:)))
                    error('DamageHomogenizationFitter:InvalidParameterVector', ...
                        'Parameter vectors must contain finite values.');
                end

                if numel(unique(p(:))) < 2
                    error('DamageHomogenizationFitter:InvalidParameterVector', ...
                        'Each parameter vector must contain at least two distinct values.');
                end

            end

        end

        function validateTensorDimensions(paramVectors, C)

            if ~isnumeric(C)
                error('DamageHomogenizationFitter:InvalidTensor', ...
                    'C must be numeric.');
            end

            for k = 1:4
                if size(C, k) ~= 2
                    error('DamageHomogenizationFitter:InvalidTensor', ...
                        'The first four dimensions of C must be 2x2x2x2.');
                end
            end

            nVar = numel(paramVectors);

            for q = 1:nVar

                expected = numel(paramVectors{q});
                actual = size(C, 4+q);

                if actual ~= expected
                    error('DamageHomogenizationFitter:TensorParameterMismatch', ...
                        ['Tensor dimension %d has size %d, but parameter %d ', ...
                         'has %d values.'], ...
                        4+q, actual, q, expected);
                end

            end

        end

        function Y = computeCholeskyTargets(C, nVar)

            maps = { ...
                [1 1 1 1], ...
                [1 1 2 2], ...
                [1 1 1 2], ...
                [2 2 1 1], ...
                [2 2 2 2], ...
                [2 2 1 2], ...
                [1 2 1 1], ...
                [1 2 2 2], ...
                [1 2 1 2]};

            nPts = 1;
            for q = 1:nVar
                nPts = nPts*size(C, 4+q);
            end

            V = zeros(nPts, 9);

            for m = 1:9
                subs = [num2cell(maps{m}), repmat({':'}, 1, nVar)];
                values = C(subs{:});
                V(:,m) = reshape(values, [], 1);
            end

            Y = zeros(nPts, 6);

            for i = 1:nPts

                Cvoigt = [ ...
                    V(i,1), V(i,2), V(i,3); ...
                    V(i,4), V(i,5), V(i,6); ...
                    V(i,7), V(i,8), V(i,9)];

                Cvoigt = 0.5*(Cvoigt + Cvoigt') + 1e-10*eye(3);
                L = chol(Cvoigt, 'lower');

                Y(i,1) = log(L(1,1));
                Y(i,2) = L(2,1);
                Y(i,3) = log(L(2,2));
                Y(i,4) = L(3,1);
                Y(i,5) = L(3,2);
                Y(i,6) = log(L(3,3));

            end

        end

        function [iTrain, iValidation] = createSplit(P, referenceDataFile)

            nPts = size(P, 1);

            if ~isempty(referenceDataFile)

                if ~exist(referenceDataFile, 'file')
                    error('DamageHomogenizationFitter:MissingReferenceData', ...
                        'Reference data file %s does not exist.', ...
                        referenceDataFile);
                end

                ref = load(referenceDataFile, 'data');

                if ~isfield(ref, 'data') || ...
                   ~isfield(ref.data, 'parameterPoints') || ...
                   ~isfield(ref.data, 'iTrain') || ...
                   ~isfield(ref.data, 'iValidation')
                    error('DamageHomogenizationFitter:InvalidReferenceData', ...
                        ['Reference file must contain data.parameterPoints, ', ...
                         'data.iTrain and data.iValidation.']);
                end

                if ~isequal(size(ref.data.parameterPoints), size(P)) || ...
                   max(abs(ref.data.parameterPoints(:) - P(:))) > 1e-12
                    error('DamageHomogenizationFitter:ReferenceGridMismatch', ...
                        'Reference data uses a different physical parameter grid.');
                end

                iTrain = ref.data.iTrain(:)';
                iValidation = ref.data.iValidation(:)';

                return

            end

            perm = randperm(nPts);
            nTrain = round(0.8*nPts);

            iTrain = perm(1:nTrain);
            iValidation = perm(nTrain+1:end);

        end

        function problem = trainNetwork(data, cfg)

            networkParams.hiddenLayers = cfg.hiddenLayers;
            networkParams.HUtype = cfg.HUtype;
            networkParams.OUtype = cfg.OUtype;
            networkParams.data = data;

            optimizerParams.type = cfg.optimizerType;
            optimizerParams.maxEpochs = cfg.maxEpochs;
            optimizerParams.learningRate = cfg.learningRate;

            if ~isempty(cfg.earlyStop)
                optimizerParams.earlyStop = cfg.earlyStop;
            end

            costParams.costType = cfg.costType;
            costParams.lambda = cfg.lambda;

            cParams.data = data;
            cParams.networkParams = networkParams;
            cParams.optimizerParams = optimizerParams;
            cParams.costParams = costParams;

            fprintf('\n');
            fprintf('Training surrogate\n');
            fprintf('Variables      : %s\n', ...
                strjoin(cfg.parameterNames, ', '));
            fprintf('Transforms     : %s\n', ...
                strjoin(cfg.transforms, ', '));
            fprintf('Feature map    : %s\n', cfg.featureMap);

            if strcmp(cfg.featureMap, 'polynomial')
                fprintf('Polynomial deg : %d\n', cfg.pol_deg);
            end

            fprintf('Activation     : %s\n', cfg.HUtype);
            fprintf('Hidden layers  : %s\n', mat2str(cfg.hiddenLayers));
            fprintf('Learning rate  : %.8g\n', cfg.learningRate);
            fprintf('Lambda         : %.8g\n', cfg.lambda);
            fprintf('Max epochs     : %d\n', cfg.maxEpochs);
            fprintf('Seed           : %d\n\n', cfg.seed);

            problem = OptimizationProblemNN(cParams);
            problem.solve();

        end

        function validateLoadedModel(paramVectors, cfg, data, savedCfg)

            nVar = numel(paramVectors);

            if ~isfield(data, 'paramVectors') || ...
               numel(data.paramVectors) ~= nVar
                error('DamageHomogenizationFitter:SavedModelMismatch', ...
                    'Saved model has incompatible parameter data.');
            end

            for q = 1:nVar

                pRequested = paramVectors{q}(:);
                pSaved = data.paramVectors{q}(:);

                if ~isequal(size(pRequested), size(pSaved)) || ...
                   max(abs(pRequested - pSaved)) > 1e-12
                    error('DamageHomogenizationFitter:SavedModelMismatch', ...
                        'Saved model uses a different grid for parameter %d.', q);
                end

            end

            if ~isequal(cfg.parameterNames, savedCfg.parameterNames)
                error('DamageHomogenizationFitter:SavedModelMismatch', ...
                    'parameterNames do not match the saved model.');
            end

            if ~isequal(cfg.transforms, savedCfg.transforms)
                error('DamageHomogenizationFitter:SavedModelMismatch', ...
                    'transforms do not match the saved model.');
            end

            if ~strcmp(cfg.featureMap, savedCfg.featureMap)
                error('DamageHomogenizationFitter:SavedModelMismatch', ...
                    'featureMap does not match the saved model.');
            end

            if strcmp(cfg.featureMap, 'polynomial') && ...
               cfg.pol_deg ~= savedCfg.pol_deg
                error('DamageHomogenizationFitter:SavedModelMismatch', ...
                    'pol_deg does not match the saved model.');
            end

        end

        function [U, dUdp] = transformParameters(P, transforms)

            nPts = size(P, 1);
            nVar = size(P, 2);

            U = zeros(nPts, nVar);
            dUdp = zeros(nPts, nVar);

            for q = 1:nVar

                p = P(:,q);
                transform = transforms{q};

                switch transform

                    case 'identity'

                        U(:,q) = p;
                        dUdp(:,q) = 1;

                    case 'log'

                        if any(p <= 0)
                            error('DamageHomogenizationFitter:LogDomain', ...
                                'Log-transformed parameters must be positive.');
                        end

                        U(:,q) = log(p);
                        dUdp(:,q) = 1./p;

                    case 'atanh'

                        if any(abs(p) >= 1)
                            error('DamageHomogenizationFitter:AtanhDomain', ...
                                'atanh-transformed parameters must satisfy |p| < 1.');
                        end

                        U(:,q) = atanh(p);
                        dUdp(:,q) = 1./(1-p.^2);

                    otherwise

                        error('DamageHomogenizationFitter:InvalidTransform', ...
                            'Unsupported transform %s.', transform);

                end

            end

        end

        function [X, dXdz] = ...
                buildFeaturesAndDerivatives(Z, featureMap, polDeg)

            nPts = size(Z, 1);
            nVar = size(Z, 2);

            switch featureMap

                case 'direct'

                    X = Z;

                    dXdz = zeros(nPts, nVar, nVar);

                    for q = 1:nVar
                        dXdz(:,q,q) = 1;
                    end

                case 'polynomial'

                    E = DamageHomogenizationFitter. ...
                        buildPolynomialExponentMatrix(nVar, polDeg);

                    nFeatures = size(E, 1);

                    X = ones(nPts, nFeatures);

                    for r = 1:nVar
                        X = X .* Z(:,r).^E(:,r)';
                    end

                    dXdz = zeros(nPts, nFeatures, nVar);

                    for q = 1:nVar

                        eq = E(:,q)';
                        active = eq > 0;

                        if ~any(active)
                            continue
                        end

                        term = ones(nPts, sum(active));

                        for r = 1:nVar

                            er = E(active,r)';

                            if r == q
                                er = er - 1;
                            end

                            term = term .* Z(:,r).^er;

                        end

                        dXdz(:,active,q) = ...
                            term .* eq(active);

                    end

                otherwise

                    error('DamageHomogenizationFitter:InvalidFeatureMap', ...
                        'Unsupported feature map %s.', featureMap);

            end

        end

        function E = buildPolynomialExponentMatrix(nVar, d)

            blocks = cell(d, 1);

            for g = 1:d
                blocks{g} = ...
                    DamageHomogenizationFitter. ...
                    generatePolynomialExponents(nVar, g);
            end

            E = vertcat(blocks{:});

        end

        function E = generatePolynomialExponents(n, g)

            if n == 1
                E = g;
                return
            end

            E = [];

            for k = 0:g

                S = DamageHomogenizationFitter. ...
                    generatePolynomialExponents(n-1, g-k);

                E = [ ...
                    E; ...
                    k*ones(size(S,1),1), S];

            end

        end

        function [fun, dfun, ddfun] = buildHandles(problem, data)

            fun = cell(2,2,2,2);
            dfun = cell(2,2,2,2);
            ddfun = cell(2,2,2,2);

            compMap = { ...
                [1 1 1 1], ...
                [2 2 2 2], ...
                [1 1 2 2], ...
                [1 2 1 2], ...
                [1 1 1 2], ...
                [2 2 1 2]};

            cache = containers.Map( ...
                'KeyType', 'char', ...
                'ValueType', 'any');

            cache('valid') = false;
            cache('P') = [];
            cache('val') = [];
            cache('dC') = [];

            for m = 1:6

                idx = compMap{m};

                i = idx(1);
                j = idx(2);
                k = idx(3);
                l = idx(4);

                f = DamageHomogenizationFitter.makeEval( ...
                    problem, data, m, cache);

                df = DamageHomogenizationFitter.makeGrad( ...
                    problem, data, m, cache);

                symmetries = { ...
                    [i,j,k,l], ...
                    [j,i,k,l], ...
                    [i,j,l,k], ...
                    [j,i,l,k], ...
                    [k,l,i,j], ...
                    [l,k,i,j], ...
                    [k,l,j,i], ...
                    [l,k,j,i]};

                for isym = 1:numel(symmetries)

                    s = symmetries{isym};

                    fun{s(1),s(2),s(3),s(4)} = f;
                    dfun{s(1),s(2),s(3),s(4)} = df;

                end

            end

        end

        function h = makeEval(problem, data, m, cache)

            h = @(varargin) ...
                DamageHomogenizationFitter.evalCompCached( ...
                problem, data, m, cache, varargin{:});

        end

        function h = makeGrad(problem, data, m, cache)

            h = @(varargin) ...
                DamageHomogenizationFitter.evalGradCached( ...
                problem, data, m, cache, varargin{:});

        end

        function val = evalCompCached( ...
                problem, data, m, cache, varargin)

            [P, outputSize] = ...
                DamageHomogenizationFitter.collectInputs( ...
                data.nVariables, varargin{:});

            DamageHomogenizationFitter.refreshCache( ...
                problem, data, P, cache);

            valAll = cache('val');

            val = reshape( ...
                valAll(:,m), ...
                outputSize);

        end

        function dJ = evalGradCached( ...
                problem, data, m, cache, varargin)

            [P, outputSize] = ...
                DamageHomogenizationFitter.collectInputs( ...
                data.nVariables, varargin{:});

            DamageHomogenizationFitter.refreshCache( ...
                problem, data, P, cache);

            dC = cache('dC');

            dJ = cell(1, data.nVariables);

            for q = 1:data.nVariables
                dJ{q} = reshape( ...
                    dC(:,m,q), ...
                    outputSize);
            end

        end

        function [P, outputSize] = collectInputs(nVar, varargin)

            if numel(varargin) ~= nVar
                error('DamageHomogenizationFitter:InputCount', ...
                    'Expected %d parameter inputs, received %d.', ...
                    nVar, numel(varargin));
            end

            values = cell(1, nVar);
            outputSize = [];

            for q = 1:nVar

                v = varargin{q};

                if isa(v, 'LagrangianFunction')
                    v = v.fValues;
                end

                if ~isnumeric(v)
                    error('DamageHomogenizationFitter:InvalidEvaluationInput', ...
                        'Surrogate inputs must be numeric or LagrangianFunction.');
                end

                if any(~isfinite(v(:)))
                    error('DamageHomogenizationFitter:InvalidEvaluationInput', ...
                        'Surrogate inputs must contain finite values.');
                end

                values{q} = v;

                if ~isscalar(v)

                    if isempty(outputSize)
                        outputSize = size(v);
                    elseif ~isequal(size(v), outputSize)
                        error('DamageHomogenizationFitter:InputSizeMismatch', ...
                            'All non-scalar surrogate inputs must have the same size.');
                    end

                end

            end

            if isempty(outputSize)
                outputSize = [1 1];
            end

            nPts = prod(outputSize);
            P = zeros(nPts, nVar);

            for q = 1:nVar

                v = values{q};

                if isscalar(v)
                    P(:,q) = v;
                else
                    P(:,q) = v(:);
                end

            end

        end

        function refreshCache(problem, data, P, cache)

            valid = cache('valid');
            oldP = cache('P');

            cacheValid = false;

            if valid && isequal(size(oldP), size(P))

                if isempty(P)
                    cacheValid = true;
                else
                    cacheValid = ...
                        max(abs(oldP(:) - P(:))) < 1e-12;
                end

            end

            if ~cacheValid

                [valAll, dC] = ...
                    DamageHomogenizationFitter.evalAllComponents( ...
                    problem, data, P);

                cache('valid') = true;
                cache('P') = P;
                cache('val') = valAll;
                cache('dC') = dC;

            end

        end

        function [valAll, dC] = ...
                evalAllComponents(problem, data, P)

            DamageHomogenizationFitter.validateEvaluationDomain( ...
                P, data);

            [U, dUdp] = ...
                DamageHomogenizationFitter.transformParameters( ...
                P, data.transforms);

            uRange = data.uMax - data.uMin;

            Z = 2*(U - data.uMin)./uRange - 1;

            dZdp = dUdp .* (2./uRange);

            [Xraw, dXdz] = ...
                DamageHomogenizationFitter.buildFeaturesAndDerivatives( ...
                Z, data.featureMap, data.pol_deg);

            Xn = (Xraw - data.muX)./data.stdX;

            nPts = size(P, 1);
            nVar = data.nVariables;
            nFeatures = data.nFeatures;

            dX = zeros(nPts, nFeatures, nVar);

            for q = 1:nVar

                dFeatureDp = ...
                    dXdz(:,:,q) .* dZdp(:,q);

                dX(:,:,q) = ...
                    dFeatureDp ./ data.stdX;

            end

            Yn = problem.computeOutputValues(Xn);

            Y = Yn .* data.stdY + data.muY;

            dYn = problem.computeDirectionalGradient(Xn, dX);

            dY = zeros(nPts, 6, nVar);

            for q = 1:nVar
                dY(:,:,q) = ...
                    dYn(:,:,q) .* data.stdY;
            end

            [valAll, dC] = ...
                DamageHomogenizationFitter. ...
                reconstructTensorAndGradient(Y, dY);

        end

        function validateEvaluationDomain(P, data)

            if size(P, 2) ~= data.nVariables
                error('DamageHomogenizationFitter:InputDimension', ...
                    'Incorrect number of surrogate variables.');
            end

            if data.allowExtrapolation
                return
            end

            tol = 1e-12;

            below = P < data.pMin - tol;
            above = P > data.pMax + tol;

            if any(below(:)) || any(above(:))
                error('DamageHomogenizationFitter:Extrapolation', ...
                    'Surrogate evaluation requested outside the training domain.');
            end

        end

        function [valAll, dC] = ...
                reconstructTensorAndGradient(Y, dY)

            nPts = size(Y, 1);
            nVar = size(dY, 3);

            L11 = exp(Y(:,1));
            L21 = Y(:,2);
            L22 = exp(Y(:,3));
            L31 = Y(:,4);
            L32 = Y(:,5);
            L33 = exp(Y(:,6));

            valAll = zeros(nPts, 6);

            valAll(:,1) = L11.^2;
            valAll(:,2) = L21.^2 + L22.^2;
            valAll(:,3) = L11.*L21;
            valAll(:,4) = L31.^2 + L32.^2 + L33.^2;
            valAll(:,5) = L11.*L31;
            valAll(:,6) = L21.*L31 + L22.*L32;

            dC = zeros(nPts, 6, nVar);

            for q = 1:nVar

                dL11 = L11 .* dY(:,1,q);
                dL21 = dY(:,2,q);
                dL22 = L22 .* dY(:,3,q);
                dL31 = dY(:,4,q);
                dL32 = dY(:,5,q);
                dL33 = L33 .* dY(:,6,q);

                dC(:,1,q) = ...
                    2*L11.*dL11;

                dC(:,2,q) = ...
                    2*L21.*dL21 + ...
                    2*L22.*dL22;

                dC(:,3,q) = ...
                    dL11.*L21 + ...
                    L11.*dL21;

                dC(:,4,q) = ...
                    2*L31.*dL31 + ...
                    2*L32.*dL32 + ...
                    2*L33.*dL33;

                dC(:,5,q) = ...
                    dL11.*L31 + ...
                    L11.*dL31;

                dC(:,6,q) = ...
                    dL21.*L31 + ...
                    L21.*dL31 + ...
                    dL22.*L32 + ...
                    L22.*dL32;

            end

        end

    end

end
