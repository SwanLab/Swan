% classdef DamageHomogenizationFitter < handle
% 
%     methods (Access = public, Static)
% 
%         function [fun, dfun, ddfun] = computePolynomial(degPoly, phi, C)
%             obj = DamageHomogenizationFitter();
%             fun = obj.computeFitting(degPoly, phi, C);
%             [dfun, ddfun] = obj.computeDerivative(fun);
%             [fun, dfun, ddfun] = obj.convertToHandle(fun, dfun, ddfun);
%         end
% 
%         function [fun, dfun, ddfun] = computeNN(b_vec, rho_vec, C, varargin)
%             saveFile     = 'HomogNN.mat';
%             forceRetrain = false;
%             if nargin >= 4 && isstruct(varargin{1}) && isfield(varargin{1}, 'retrain')
%                 forceRetrain = varargin{1}.retrain;
%             end
%             if ~forceRetrain && exist(saveFile, 'file')
%                 fprintf('Carregando rede guardada de %s...\n', saveFile);
%                 loaded  = load(saveFile, 'problem', 'data');
%                 problem = loaded.problem;
%                 data    = loaded.data;
%             else
%                 data    = DamageHomogenizationFitter.prepareData(b_vec, rho_vec, C);
%                 problem = DamageHomogenizationFitter.trainNetwork(data, varargin{:});
%                 save(saveFile, 'problem', 'data');
%                 fprintf('Rede guardada em %s\n', saveFile);
%             end
%             [fun, dfun, ddfun] = DamageHomogenizationFitter.buildHandles(problem, data);
%         end
% 
%     end
% 
%     methods (Access = private)
% 
%         function fun = computeFitting(~, degPoly, phi, C)
%             phi   = reshape(phi, length(phi), []);
%             nStre = size(C, 1);
%             fun   = cell(2,2,2,2);
%             for i = 1:nStre
%                 for j = 1:nStre
%                     for k = 1:nStre
%                         for l = 1:nStre
%                             coeffs       = polyfit(phi, squeeze(C(i,j,k,l,:)), degPoly);
%                             fun{i,j,k,l} = poly2sym(coeffs);
%                             if isempty(symvar(fun{i,j,k,l}))
%                                 syms x
%                                 fun{i,j,k,l} = 1e-20.*x.^9;
%                             end
%                         end
%                     end
%                 end
%             end
%         end
% 
%         function [dfun, ddfun] = computeDerivative(~, fun)
%             nStre = size(fun, 1);
%             dfun  = cell(2,2,2,2);
%             ddfun = cell(2,2,2,2);
%             for i = 1:nStre
%                 for j = 1:nStre
%                     for k = 1:nStre
%                         for l = 1:nStre
%                             dfun{i,j,k,l}  = diff(fun{i,j,k,l});
%                             ddfun{i,j,k,l} = diff(dfun{i,j,k,l});
%                         end
%                     end
%                 end
%             end
%         end
% 
%         function [fun, dfun, ddfun] = convertToHandle(~, fun, dfun, ddfun)
%             nStre = size(fun, 1);
%             for i = 1:nStre
%                 for j = 1:nStre
%                     for k = 1:nStre
%                         for l = 1:nStre
%                             fun{i,j,k,l}   = matlabFunction(fun{i,j,k,l});
%                             dfun{i,j,k,l}  = matlabFunction(dfun{i,j,k,l});
%                             ddfun{i,j,k,l} = matlabFunction(ddfun{i,j,k,l});
%                         end
%                     end
%                 end
%             end
%         end
% 
%     end
% 
%     methods (Access = private, Static)
% 
%         function data = prepareData(b_vec, rho_vec, C)
%             nB   = length(b_vec);
%             nRho = length(rho_vec);
%             X = zeros(nB*nRho,2);
% 
%             row = 1;
%             for iRho = 1:nRho
%                 for iB = 1:nB
%                     X(row,1) = b_vec(iB);
%                     X(row,2) = rho_vec(iRho);
%                     row = row + 1;
%                 end
%             end
% 
%             % =============================================
%             % Normalizar entradas para [-1, 1] (geral)
%             % =============================================
%             b_min = min(b_vec);
%             b_max = max(b_vec);
%             rho_min = min(rho_vec);
%             rho_max = max(rho_vec);
% 
%             X_norm = X;
%             X_norm(:,1) = 2 * (X(:,1) - b_min) / (b_max - b_min) - 1;
%             X_norm(:,2) = 2 * (X(:,2) - rho_min) / (rho_max - rho_min) - 1;
%             % =============================================
% 
%             compMap = {[1,1,1,1],[2,2,2,2],[1,1,2,2],[1,2,1,2],[1,1,1,2],[2,2,1,2]};
%             nComp = length(compMap);
%             nPts  = size(X, 1);
%             Y = zeros(nPts, nComp);
%             for iRho = 1:nRho
%                 for iB = 1:nB
%                     row = (iRho-1)*nB + iB;
%                     for m = 1:nComp
%                         idx = compMap{m};
%                         Y(row, m) = C(idx(1),idx(2),idx(3),idx(4),iRho,iB);
%                     end
%                 end
%             end
% 
%             muX  = mean(X_norm, 1);
%             stdX = std(X_norm, 0, 1);
%             muY  = mean(Y, 1);
%             stdY = std(Y, 0, 1);
%             stdX(stdX == 0) = 1;
%             stdY(stdY == 0) = 1;
% 
%             Xn = (X_norm - muX) ./ stdX;
%             Yn = (Y - muY) ./ stdY;
% 
%             perm   = randperm(nPts);
%             nTrain = round(0.8 * nPts);
%             iTrain = perm(1:nTrain);
%             iTest  = perm(nTrain+1:end);
% 
%             data.Xtrain    = Xn(iTrain, :);
%             data.Ytrain    = Yn(iTrain, :);
%             data.Xtest     = Xn(iTest,  :);
%             data.Ytest     = Yn(iTest,  :);
%             data.nFeatures = 2;
%             data.nLabels   = nComp;
%             data.muX       = muX;
%             data.stdX      = stdX;
%             data.muY       = muY;
%             data.stdY      = stdY;
%             data.compMap   = compMap;
% 
%             % Guardar parâmetros de normalização
%             data.b_min = b_min;
%             data.b_max = b_max;
%             data.rho_min = rho_min;
%             data.rho_max = rho_max;
%         end
% 
%         function problem = trainNetwork(data, varargin)
%             hiddenLayers = [64, 64, 32];
%             maxEpochs    = 4000;
%             learningRate = 1e-3;
%             if nargin >= 2 && isstruct(varargin{1})
%                 p = varargin{1};
%                 if isfield(p,'hiddenLayers'), hiddenLayers = p.hiddenLayers; end
%                 if isfield(p,'maxEpochs'),    maxEpochs    = p.maxEpochs;    end
%                 if isfield(p,'learningRate'), learningRate = p.learningRate; end
%             end
%             networkParams.hiddenLayers = hiddenLayers;
%             networkParams.HUtype       = 'tanh';
%             networkParams.OUtype       = 'linear';
%             networkParams.data         = data;
%             optimizerParams.type         = 'SGD';
%             optimizerParams.maxEpochs    = maxEpochs;
%             optimizerParams.learningRate = learningRate;
%             costParams.costType = 'L2';
%             costParams.lambda   = 1e-3;
%             cParams.data            = data;
%             cParams.networkParams   = networkParams;
%             cParams.optimizerParams = optimizerParams;
%             cParams.costParams      = costParams;
%             problem = OptimizationProblemNN(cParams);
%             problem.solve();
%         end
% 
%         function [fun, dfun, ddfun] = buildHandles(problem, data)
%             compMap = data.compMap;
%             nComp   = length(compMap);
%             fun     = cell(2,2,2,2);
%             dfun    = cell(2,2,2,2);
%             ddfun   = cell(2,2,2,2);
% 
%             for m = 1:nComp
%                 idx = compMap{m};
%                 i=idx(1); j=idx(2); k=idx(3); l=idx(4);
% 
%                 f  = DamageHomogenizationFitter.makeEval(problem, data, m);
%                 df = DamageHomogenizationFitter.makeGrad(problem, data, m);
%                 dd = DamageHomogenizationFitter.makeHess(problem, data, m);
% 
% 
%                 for si = {[i,j,k,l],[j,i,k,l],[i,j,l,k],[j,i,l,k], ...
%                            [k,l,i,j],[l,k,i,j],[k,l,j,i],[l,k,j,i]}
%                     s = si{1};
%                     fun{s(1),s(2),s(3),s(4)}   = f;
%                     dfun{s(1),s(2),s(3),s(4)}  = df;
%                     ddfun{s(1),s(2),s(3),s(4)} = dd;
%                 end
%             end
%         end
% 
%         function h = makeEval(problem, data, m)
%             h = @(b, rho) DamageHomogenizationFitter.evalComp(problem, data, b, rho, m);
%         end
% 
%         function val = evalComp(problem, data, b, rho, m)
%             if isa(b,   'LagrangianFunction'), b   = b.fValues; end
%             if isa(rho, 'LagrangianFunction'), rho = rho.fValues; end
%             origSize = size(b);
%             b   = b(:);
%             rho = rho(:);
% 
%             % Normalizar entrada
%             b_norm = 2 * (b - data.b_min) / (data.b_max - data.b_min) - 1;
%             rho_norm = 2 * (rho - data.rho_min) / (data.rho_max - data.rho_min) - 1;
% 
%             Xn = ([b_norm, rho_norm] - data.muX) ./ data.stdX;
%             Yn = problem.computeOutputValues(Xn);
%             Y   = Yn .* data.stdY + data.muY;
%             val = reshape(Y(:, m), origSize);
%         end
% 
%         function h = makeGrad(problem, data, m)
%             h = @(b, rho) DamageHomogenizationFitter.evalGrad(problem, data, b, rho, m);
%         end
% 
%        function dJ = evalGrad(problem, data, b, rho, m)
%             if isa(b,   'LagrangianFunction'), b   = b.fValues; end
%             if isa(rho, 'LagrangianFunction'), rho = rho.fValues; end
% 
%             origSize = size(b);
%             b   = b(:);
%             rho = rho(:);
% 
%             % Normalizar entrada
%             b_norm = 2 * (b - data.b_min) / (data.b_max - data.b_min) - 1;
%             rho_norm = 2 * (rho - data.rho_min) / (data.rho_max - data.rho_min) - 1;
% 
%             Xn = ([b_norm, rho_norm] - data.muX) ./ data.stdX;
% 
%             % J: (nPts, nLabels, nFeatures)
%             J = problem.computeGradient(Xn);
% 
%             scale = data.stdY(m);
% 
%             % Gradiente em relação à entrada normalizada
%             dJ_db_norm = J(:, m, 1) * (scale / data.stdX(1));
%             dJ_drho_norm = J(:, m, 2) * (scale / data.stdX(2));
% 
%             % Chain rule para desnormalizar
%             % db_norm/db = 2/(b_max - b_min)
%             % drho_norm/drho = 2/(rho_max - rho_min)
%             dJ = cell(1, 2);
%             dJ{1} = reshape(dJ_db_norm * 2 / (data.b_max - data.b_min), origSize);
%             dJ{2} = reshape(dJ_drho_norm * 2 / (data.rho_max - data.rho_min), origSize);
%         end
%         function h = makeHess(problem, data, m)
%             h = @(b, rho) DamageHomogenizationFitter.evalHess(problem, data, b, rho, m);
%         end
% 
%         function ddJ = evalHess(problem, data, b, rho, m)
%             eps = 1e-5;
%             dJ_pb  = DamageHomogenizationFitter.evalGrad(problem, data, b+eps, rho,     m);
%             dJ_mb  = DamageHomogenizationFitter.evalGrad(problem, data, b-eps, rho,     m);
%             dJ_pr  = DamageHomogenizationFitter.evalGrad(problem, data, b,     rho+eps, m);
%             dJ_mr  = DamageHomogenizationFitter.evalGrad(problem, data, b,     rho-eps, m);
%             ddJ = cell(2,2);
%             ddJ{1,1} = (dJ_pb{1} - dJ_mb{1}) / (2*eps);
%             ddJ{2,2} = (dJ_pr{2} - dJ_mr{2}) / (2*eps);
%             ddJ{1,2} = (dJ_pr{1} - dJ_mr{1}) / (2*eps);
%             ddJ{2,1} = ddJ{1,2};
%         end
% 
%     end
% 
% end
% 
classdef DamageHomogenizationFitter < handle

    methods (Access = public, Static)

        function [fun, dfun, ddfun] = computePolynomial(degPoly, phi, C)
            obj = DamageHomogenizationFitter();
            fun = obj.computeFitting(degPoly, phi, C);
            [dfun, ddfun] = obj.computeDerivative(fun);
            [fun, dfun, ddfun] = obj.convertToHandle(fun, dfun, ddfun);
        end

        function [fun, dfun, ddfun] = computeNN(b_vec, rho_vec, C, varargin)

            saveFile     = 'HomogNN.mat';
            forceRetrain = false;
            pol_deg      = 1;

            if nargin >= 4 && isstruct(varargin{1})
                p = varargin{1};
                if isfield(p, 'retrain'),  forceRetrain = p.retrain;  end
                if isfield(p, 'pol_deg'),  pol_deg      = p.pol_deg;  end
            end

            if ~forceRetrain && exist(saveFile, 'file')
                fprintf('Carregando rede guardada de %s...\n', saveFile);
                loaded  = load(saveFile, 'problem', 'data');
                problem = loaded.problem;
                data    = loaded.data;

                precisaRetreinar = ~isfield(data, 'pol_deg') || data.pol_deg ~= pol_deg;
                if precisaRetreinar
                    if ~isfield(data, 'pol_deg')
                        warning('%s versao antiga sem pol_deg. Retreinando com pol_deg=%d...', saveFile, pol_deg);
                    else
                        warning('%s treinado com pol_deg=%d, solicitado pol_deg=%d. Retreinando...', saveFile, data.pol_deg, pol_deg);
                    end
                    data    = DamageHomogenizationFitter.prepareData(b_vec, rho_vec, C, pol_deg);
                    problem = DamageHomogenizationFitter.trainNetwork(data, varargin{:});
                    history = problem.getHistory();
                    save(saveFile, 'problem', 'data', 'history');
                    save('NNhistory.mat', 'history');
                end
            else
                data    = DamageHomogenizationFitter.prepareData(b_vec, rho_vec, C, pol_deg);
                problem = DamageHomogenizationFitter.trainNetwork(data, varargin{:});
                history = problem.getHistory();
                save(saveFile, 'problem', 'data', 'history');
                save('NNhistory.mat', 'history');
                fprintf('Rede guardada em %s\n', saveFile);
            end

            [fun, dfun, ddfun] = DamageHomogenizationFitter.buildHandles(problem, data);
        end

    end

    methods (Access = private)

        function fun = computeFitting(~, degPoly, phi, C)
            phi   = reshape(phi, length(phi), []);
            nStre = size(C, 1);
            fun   = cell(2,2,2,2);
            for i = 1:nStre
                for j = 1:nStre
                    for k = 1:nStre
                        for l = 1:nStre
                            coeffs       = polyfit(phi, squeeze(C(i,j,k,l,:)), degPoly);
                            fun{i,j,k,l} = poly2sym(coeffs);
                            if isempty(symvar(fun{i,j,k,l}))
                                syms x
                                fun{i,j,k,l} = 1e-20.*x.^9;
                            end
                        end
                    end
                end
            end
        end

        function [dfun, ddfun] = computeDerivative(~, fun)
            nStre = size(fun, 1);
            dfun  = cell(2,2,2,2);
            ddfun = cell(2,2,2,2);
            for i = 1:nStre
                for j = 1:nStre
                    for k = 1:nStre
                        for l = 1:nStre
                            dfun{i,j,k,l}  = diff(fun{i,j,k,l});
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
                            fun{i,j,k,l}   = matlabFunction(fun{i,j,k,l});
                            dfun{i,j,k,l}  = matlabFunction(dfun{i,j,k,l});
                            ddfun{i,j,k,l} = matlabFunction(ddfun{i,j,k,l});
                        end
                    end
                end
            end
        end

    end

    methods (Access = private, Static)

        function data = prepareData(b_vec, rho_vec, C, pol_deg)
            if nargin < 4, pol_deg = 1; end

            nB   = length(b_vec);
            nRho = length(rho_vec);
            X    = zeros(nB*nRho, 2);
            row  = 1;
            for iRho = 1:nRho
                for iB = 1:nB
                    X(row,1) = b_vec(iB);
                    X(row,2) = rho_vec(iRho);
                    row = row + 1;
                end
            end

            b_min   = min(b_vec);   b_max   = max(b_vec);
            rho_min = min(rho_vec); rho_max = max(rho_vec);

            X_norm      = X;
            X_norm(:,1) = 2*(X(:,1)-b_min)/(b_max-b_min) - 1;
            X_norm(:,2) = 2*(X(:,2)-rho_min)/(rho_max-rho_min) - 1;

            X_poly = DamageHomogenizationFitter.buildPolynomialFeatures(X_norm, pol_deg);

            nPts = size(X,1);
            Y    = zeros(nPts, 6);
            for iRho = 1:nRho
                for iB = 1:nB
                    row = (iRho-1)*nB + iB;
                    Cvoigt = [
                        C(1,1,1,1,iRho,iB), C(1,1,2,2,iRho,iB), C(1,1,1,2,iRho,iB);
                        C(2,2,1,1,iRho,iB), C(2,2,2,2,iRho,iB), C(2,2,1,2,iRho,iB);
                        C(1,2,1,1,iRho,iB), C(1,2,2,2,iRho,iB), C(1,2,1,2,iRho,iB)];
                    Cvoigt = 0.5*(Cvoigt+Cvoigt') + 1e-10*eye(3);
                    L = chol(Cvoigt,'lower');
                    Y(row,1) = log(L(1,1));
                    Y(row,2) = L(2,1);
                    Y(row,3) = log(L(2,2));
                    Y(row,4) = L(3,1);
                    Y(row,5) = L(3,2);
                    Y(row,6) = log(L(3,3));
                end
            end

            muX  = mean(X_poly,1); stdX = std(X_poly,0,1);
            muY  = mean(Y,1);      stdY = std(Y,0,1);
            stdX(stdX==0) = 1;     stdY(stdY==0) = 1;

            Xn = (X_poly-muX)./stdX;
            Yn = (Y-muY)./stdY;

            perm   = randperm(nPts);
            nTrain = round(0.8*nPts);
            iTrain = perm(1:nTrain);
            iTest  = perm(nTrain+1:end);

            data.Xtrain    = Xn(iTrain,:);  data.Ytrain  = Yn(iTrain,:);
            data.Xtest     = Xn(iTest,:);   data.Ytest   = Yn(iTest,:);
            data.nFeatures = size(Xn,2);    data.nLabels = 6;
            data.muX = muX;  data.stdX = stdX;
            data.muY = muY;  data.stdY = stdY;
            data.compMap   = {1,2,3,4,5,6};
            data.b_min     = b_min;   data.b_max   = b_max;
            data.rho_min   = rho_min; data.rho_max = rho_max;
            data.pol_deg   = pol_deg;
        end

        function Xpoly = buildPolynomialFeatures(X, d)
            [N, n] = size(X);
            Xpoly  = [];
            for g = 1:d
                E = DamageHomogenizationFitter.generatePolynomialExponents(n, g);
                for k = 1:size(E,1)
                    c = ones(N,1);
                    for j = 1:n, c = c .* X(:,j).^E(k,j); end
                    Xpoly = [Xpoly, c]; %#ok<AGROW>
                end
            end
        end

        function E = generatePolynomialExponents(n, g)
            if n == 1, E = g; return; end
            E = [];
            for k = 0:g
                S = DamageHomogenizationFitter.generatePolynomialExponents(n-1, g-k);
                E = [E; k*ones(size(S,1),1), S]; %#ok<AGROW>
            end
        end

        function [dXdb, dXdrho] = buildPolynomialFeatureGradients(X, d)
            [N, n] = size(X);
            if n ~= 2, error('Assume 2 variaveis.'); end
            x1 = X(:,1); x2 = X(:,2);
            dXdb = []; dXdrho = [];
            for g = 1:d
                E = DamageHomogenizationFitter.generatePolynomialExponents(n, g);
                for k = 1:size(E,1)
                    e1 = E(k,1); e2 = E(k,2);
                    db = zeros(N,1);
                    dr = zeros(N,1);
                    if e1 > 0, db = e1 .* x1.^(e1-1) .* x2.^e2;  end
                    if e2 > 0, dr = x1.^e1 .* e2 .* x2.^(e2-1); end
                    dXdb   = [dXdb,   db]; %#ok<AGROW>
                    dXdrho = [dXdrho, dr]; %#ok<AGROW>
                end
            end
        end

        function problem = trainNetwork(data, varargin)
                
            hiddenLayers = [50 100 200 100 50 30];
            maxEpochs    = 500000;
            learningRate = 0.015;

            if nargin >= 2 && isstruct(varargin{1})
                p = varargin{1};
                if isfield(p,'hiddenLayers'), hiddenLayers = p.hiddenLayers; end
                if isfield(p,'maxEpochs'),    maxEpochs    = p.maxEpochs;    end
                if isfield(p,'learningRate'), learningRate = p.learningRate; end
            end
            networkParams.hiddenLayers = hiddenLayers;
            networkParams.HUtype       = 'ReLU';
            networkParams.OUtype       = 'linear';
            networkParams.data         = data;
            optimizerParams.type         = 'SGD';
            optimizerParams.maxEpochs    = maxEpochs;
            optimizerParams.learningRate = learningRate;
            costParams.costType = 'L2';
            costParams.lambda   = 0;
            cParams.data            = data;
            cParams.networkParams   = networkParams;
            cParams.optimizerParams = optimizerParams;
            cParams.costParams      = costParams;
            problem = OptimizationProblemNN(cParams);
            problem.solve();
        end

        % -----------------------------------------------------------------
        % buildHandles: cache via containers.Map (handle semantics),
        % garantindo que todas as anonymous functions partilham o mesmo
        % objeto de cache por referencia -- nao por valor como as celulas.
        % ddfun = vazio pois OptimizerNullSpace nao usa o Hessiano.
        % -----------------------------------------------------------------
        function [fun, dfun, ddfun] = buildHandles(problem, data)
            fun   = cell(2,2,2,2);
            dfun  = cell(2,2,2,2);
            ddfun = cell(2,2,2,2);   % vazio -- nao usado pelo otimizador

            compMap = {[1,1,1,1],[2,2,2,2],[1,1,2,2],[1,2,1,2],[1,1,1,2],[2,2,1,2]};

            % containers.Map tem semantica de handle em MATLAB:
            % e capturado POR REFERENCIA pelas anonymous functions,
            % entao todas as funcoes partilham o mesmo cache.
            cache = containers.Map('KeyType','char','ValueType','any');
            cache('valid') = false;
            cache('b')     = [];
            cache('rho')   = [];
            cache('val')   = [];
            cache('dCdb')  = [];
            cache('dCdrho')= [];

            for m = 1:6
                idx = compMap{m};
                i = idx(1); j = idx(2); k = idx(3); l = idx(4);

                f  = DamageHomogenizationFitter.makeEval(problem, data, m, cache);
                df = DamageHomogenizationFitter.makeGrad(problem, data, m, cache);

                for si = {[i,j,k,l],[j,i,k,l],[i,j,l,k],[j,i,l,k], ...
                           [k,l,i,j],[l,k,i,j],[k,l,j,i],[l,k,j,i]}
                    s = si{1};
                    fun{s(1),s(2),s(3),s(4)}  = f;
                    dfun{s(1),s(2),s(3),s(4)} = df;
                end
            end
        end

        function h = makeEval(problem, data, m, cache)
            h = @(b,rho) DamageHomogenizationFitter.evalCompCached(problem, data, b, rho, m, cache);
        end

        function h = makeGrad(problem, data, m, cache)
            h = @(b,rho) DamageHomogenizationFitter.evalGradCached(problem, data, b, rho, m, cache);
        end

        function val = evalCompCached(problem, data, b, rho, m, cache)
            if isa(b,   'LagrangianFunction'), b   = b.fValues; end
            if isa(rho, 'LagrangianFunction'), rho = rho.fValues; end
            origSize = size(b);
            bv   = b(:);
            rhov = rho(:);
            DamageHomogenizationFitter.refreshCache(problem, data, bv, rhov, cache);
            val_all = cache('val');
            val = reshape(val_all(:,m), origSize);
        end

        function dJ = evalGradCached(problem, data, b, rho, m, cache)
            if isa(b,   'LagrangianFunction'), b   = b.fValues; end
            if isa(rho, 'LagrangianFunction'), rho = rho.fValues; end
            origSize = size(b);
            bv   = b(:);
            rhov = rho(:);
            DamageHomogenizationFitter.refreshCache(problem, data, bv, rhov, cache);
            dJ = cell(1,2);
            dCdb_c   = cache('dCdb');
            dCdrho_c = cache('dCdrho');
            dJ{1} = reshape(dCdb_c(:,m),   origSize);
            dJ{2} = reshape(dCdrho_c(:,m), origSize);
        end

        function refreshCache(problem, data, bv, rhov, cache)
            tol    = 1e-12;
            valid  = cache('valid');
            bOld   = cache('b');
            rhoOld = cache('rho');
            cacheValid = valid                                && ...
                         isequal(size(bOld), size(bv))       && ...
                         max(abs(bOld   - bv))   < tol       && ...
                         max(abs(rhoOld - rhov)) < tol;

            if ~cacheValid
                [dCdb_all, dCdrho_all, val_all] = ...
                    DamageHomogenizationFitter.evalAllComponents(problem, data, bv, rhov);
                cache('valid')  = true;
                cache('b')      = bv;
                cache('rho')    = rhov;
                cache('val')    = val_all;
                cache('dCdb')   = dCdb_all;
                cache('dCdrho') = dCdrho_all;
            end
        end

        function [dCdb_all, dCdrho_all, val_all] = evalAllComponents(problem, data, b, rho)
            b_norm   = 2*(b - data.b_min)/(data.b_max - data.b_min) - 1;
            rho_norm = 2*(rho - data.rho_min)/(data.rho_max - data.rho_min) - 1;

            X_base = [b_norm, rho_norm];
            X_poly = DamageHomogenizationFitter.buildPolynomialFeatures(X_base, data.pol_deg);
            Xn     = (X_poly - data.muX) ./ data.stdX;

            Yn = problem.computeOutputValues(Xn);
            J  = problem.computeGradient(Xn);
            Y  = Yn .* data.stdY + data.muY;

            L11 = exp(Y(:,1)); L21 = Y(:,2); L22 = exp(Y(:,3));
            L31 = Y(:,4);      L32 = Y(:,5); L33 = exp(Y(:,6));

            nPts = size(b,1);
            val_all = zeros(nPts, 6);
            val_all(:,1) = L11.^2;
            val_all(:,2) = L21.^2 + L22.^2;
            val_all(:,3) = L11.*L21;
            val_all(:,4) = L31.^2 + L32.^2 + L33.^2;
            val_all(:,5) = L11.*L31;
            val_all(:,6) = L21.*L31 + L22.*L32;

            [dXpoly_dbn, dXpoly_drn] = ...
                DamageHomogenizationFitter.buildPolynomialFeatureGradients(X_base, data.pol_deg);

            nFeat = size(X_poly, 2);
            dY_db   = zeros(nPts, 6);
            dY_drho = zeros(nPts, 6);

            for mOut = 1:6
                for q = 1:nFeat
                    Jq = J(:,mOut,q) ./ data.stdX(q);
                    dY_db(:,mOut)   = dY_db(:,mOut)   + Jq .* dXpoly_dbn(:,q);
                    dY_drho(:,mOut) = dY_drho(:,mOut) + Jq .* dXpoly_drn(:,q);
                end
                dY_db(:,mOut)   = dY_db(:,mOut)   .* data.stdY(mOut) .* (2/(data.b_max   - data.b_min));
                dY_drho(:,mOut) = dY_drho(:,mOut) .* data.stdY(mOut) .* (2/(data.rho_max - data.rho_min));
            end

            dL11_db = L11.*dY_db(:,1);   dL11_drho = L11.*dY_drho(:,1);
            dL21_db = dY_db(:,2);         dL21_drho = dY_drho(:,2);
            dL22_db = L22.*dY_db(:,3);   dL22_drho = L22.*dY_drho(:,3);
            dL31_db = dY_db(:,4);         dL31_drho = dY_drho(:,4);
            dL32_db = dY_db(:,5);         dL32_drho = dY_drho(:,5);
            dL33_db = L33.*dY_db(:,6);   dL33_drho = L33.*dY_drho(:,6);

            dCdb_all   = zeros(nPts, 6);
            dCdrho_all = zeros(nPts, 6);

            dCdb_all(:,1) = 2*L11.*dL11_db;
            dCdb_all(:,2) = 2*L21.*dL21_db + 2*L22.*dL22_db;
            dCdb_all(:,3) = dL11_db.*L21 + L11.*dL21_db;
            dCdb_all(:,4) = 2*L31.*dL31_db + 2*L32.*dL32_db + 2*L33.*dL33_db;
            dCdb_all(:,5) = dL11_db.*L31 + L11.*dL31_db;
            dCdb_all(:,6) = dL21_db.*L31 + L21.*dL31_db + dL22_db.*L32 + L22.*dL32_db;

            dCdrho_all(:,1) = 2*L11.*dL11_drho;
            dCdrho_all(:,2) = 2*L21.*dL21_drho + 2*L22.*dL22_drho;
            dCdrho_all(:,3) = dL11_drho.*L21 + L11.*dL21_drho;
            dCdrho_all(:,4) = 2*L31.*dL31_drho + 2*L32.*dL32_drho + 2*L33.*dL33_drho;
            dCdrho_all(:,5) = dL11_drho.*L31 + L11.*dL31_drho;
            dCdrho_all(:,6) = dL21_drho.*L31 + L21.*dL31_drho + dL22_drho.*L32 + L22.*dL32_drho;
        end

    end

end