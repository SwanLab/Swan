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
            if nargin >= 4 && isstruct(varargin{1}) && isfield(varargin{1}, 'retrain')
                forceRetrain = varargin{1}.retrain;
            end
            if ~forceRetrain && exist(saveFile, 'file')
                fprintf('Carregando rede guardada de %s...\n', saveFile);
                loaded  = load(saveFile, 'problem', 'data');
                problem = loaded.problem;
                data    = loaded.data;
            else
                data    = DamageHomogenizationFitter.prepareData(b_vec, rho_vec, C);
                problem = DamageHomogenizationFitter.trainNetwork(data, varargin{:});
                save(saveFile, 'problem', 'data');
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

        function data = prepareData(b_vec, rho_vec, C)
            nB   = length(b_vec);
            nRho = length(rho_vec);
            X = zeros(nB*nRho,2);

            row = 1;

            for iRho = 1:nRho
                for iB = 1:nB

                    X(row,1) = b_vec(iB);
                    X(row,2) = rho_vec(iRho);

                    row = row + 1;

                end
            end
            compMap = {[1,1,1,1],[2,2,2,2],[1,1,2,2],[1,2,1,2],[1,1,1,2],[2,2,1,2]};
            nComp = length(compMap);
            nPts  = size(X, 1);
            Y = zeros(nPts, nComp);
            for iRho = 1:nRho
                for iB = 1:nB
                    row = (iRho-1)*nB + iB;
                    for m = 1:nComp
                        idx = compMap{m};
                        Y(row, m) = C(idx(1),idx(2),idx(3),idx(4),iRho,iB);
                    end
                end
            end
            muX  = mean(X, 1);  stdX = std(X, 0, 1);
            muY  = mean(Y, 1);  stdY = std(Y, 0, 1);
            stdX(stdX == 0) = 1;
            stdY(stdY == 0) = 1;
            Xn = (X - muX) ./ stdX;
            Yn = (Y - muY) ./ stdY;
            perm   = randperm(nPts);
            nTrain = round(0.8 * nPts);
            iTrain = perm(1:nTrain);
            iTest  = perm(nTrain+1:end);
            data.Xtrain    = Xn(iTrain, :);
            data.Ytrain    = Yn(iTrain, :);
            data.Xtest     = Xn(iTest,  :);
            data.Ytest     = Yn(iTest,  :);
            data.nFeatures = 2;
            data.nLabels   = nComp;
            data.muX       = muX;   data.stdX = stdX;
            data.muY       = muY;   data.stdY = stdY;
            data.compMap   = compMap;
        end

        function problem = trainNetwork(data, varargin)
            hiddenLayers = [64, 64, 32];
            maxEpochs    = 4000;
            learningRate = 1e-3;
            if nargin >= 2 && isstruct(varargin{1})
                p = varargin{1};
                if isfield(p,'hiddenLayers'), hiddenLayers = p.hiddenLayers; end
                if isfield(p,'maxEpochs'),    maxEpochs    = p.maxEpochs;    end
                if isfield(p,'learningRate'), learningRate = p.learningRate; end
            end
            networkParams.hiddenLayers = hiddenLayers;
            networkParams.HUtype       = 'tanh';
            networkParams.OUtype       = 'linear';
            networkParams.data         = data;
            optimizerParams.type         = 'SGD';
            optimizerParams.maxEpochs    = maxEpochs;
            optimizerParams.learningRate = learningRate;
            costParams.costType = 'L2';
            costParams.lambda   = 1e-4;
            cParams.data            = data;
            cParams.networkParams   = networkParams;
            cParams.optimizerParams = optimizerParams;
            cParams.costParams      = costParams;
            problem = OptimizationProblemNN(cParams);
            problem.solve();
        end

        function [fun, dfun, ddfun] = buildHandles(problem, data)
            compMap = data.compMap;
            nComp   = length(compMap);
            fun     = cell(2,2,2,2);
            dfun    = cell(2,2,2,2);
            ddfun   = cell(2,2,2,2);

            for m = 1:nComp
                idx = compMap{m};
                i=idx(1); j=idx(2); k=idx(3); l=idx(4);

                f  = DamageHomogenizationFitter.makeEval(problem, data, m);
                df = DamageHomogenizationFitter.makeGrad(problem, data, m);
                dd = DamageHomogenizationFitter.makeHess(problem, data, m);

                
                for si = {[i,j,k,l],[j,i,k,l],[i,j,l,k],[j,i,l,k], ...
                           [k,l,i,j],[l,k,i,j],[k,l,j,i],[l,k,j,i]}
                    s = si{1};
                    fun{s(1),s(2),s(3),s(4)}   = f;
                    dfun{s(1),s(2),s(3),s(4)}  = df;
                    ddfun{s(1),s(2),s(3),s(4)} = dd;
                end
            end
        end

        function h = makeEval(problem, data, m)
            h = @(b, rho) DamageHomogenizationFitter.evalComp(problem, data, b, rho, m);
        end

        function val = evalComp(problem, data, b, rho, m)
            if isa(b,   'LagrangianFunction'), b   = b.fValues; end
            if isa(rho, 'LagrangianFunction'), rho = rho.fValues; end
            origSize = size(b);
            b   = b(:);
            rho = rho(:);
            Xn  = ([b, rho] - data.muX) ./ data.stdX;
            Yn  = problem.computeOutputValues(Xn);
            Y   = Yn .* data.stdY + data.muY;
            val = reshape(Y(:, m), origSize);
        end

        function h = makeGrad(problem, data, m)
            h = @(b, rho) DamageHomogenizationFitter.evalGrad(problem, data, b, rho, m);
        end

       function dJ = evalGrad(problem, data, b, rho, m)
            if isa(b,   'LagrangianFunction'), b   = b.fValues; end
            if isa(rho, 'LagrangianFunction'), rho = rho.fValues; end

            origSize = size(b);
            b   = b(:);
            rho = rho(:);

            Xn = ([b, rho] - data.muX) ./ data.stdX;

            % J: (nPts, nLabels, nFeatures)
            J = problem.computeGradient(Xn);

            scale = data.stdY(m);
            dJ = cell(1, 2);
            dJ{1} = reshape(J(:, m, 1) * (scale / data.stdX(1)), origSize);
            dJ{2} = reshape(J(:, m, 2) * (scale / data.stdX(2)), origSize);
        end
        function h = makeHess(problem, data, m)
            h = @(b, rho) DamageHomogenizationFitter.evalHess(problem, data, b, rho, m);
        end

        function ddJ = evalHess(problem, data, b, rho, m)
            eps = 1e-5;
            dJ_pb  = DamageHomogenizationFitter.evalGrad(problem, data, b+eps, rho,     m);
            dJ_mb  = DamageHomogenizationFitter.evalGrad(problem, data, b-eps, rho,     m);
            dJ_pr  = DamageHomogenizationFitter.evalGrad(problem, data, b,     rho+eps, m);
            dJ_mr  = DamageHomogenizationFitter.evalGrad(problem, data, b,     rho-eps, m);
            ddJ = cell(2,2);
            ddJ{1,1} = (dJ_pb{1} - dJ_mb{1}) / (2*eps);
            ddJ{2,2} = (dJ_pr{2} - dJ_mr{2}) / (2*eps);
            ddJ{1,2} = (dJ_pr{1} - dJ_mr{1}) / (2*eps);
            ddJ{2,1} = ddJ{1,2};
        end

    end

end

