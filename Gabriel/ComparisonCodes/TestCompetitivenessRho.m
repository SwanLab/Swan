%% TestCompetitivenessRho.m
% Analisa a competitividade de rho no vademecum C(b,rho).
% Compara componentes normalizadas com a resposta linear y = rho/rhoMax.
%
% O script procura:
%   1) Homogenizationtwovariables15.mat na pasta atual
%   2) TOVademecum/Interpolation/Homogenizationtwovariables15.mat
%
% Componentes analisadas:
%   C1111, C2222, C1122 e C1212
%
% IMPORTANTE:
% O vademecum usado no tutorial vai ate rhoMax = 0.97; portanto, a reta
% linear correta na normalizacao por C(b,rhoMax) e rho/rhoMax, nao rho puro.

clearvars;
close all;
clc;

%% Configuracao
fileBase = 'Homogenizationtwovariables15.mat';

candidateFiles = {
    fileBase
    fullfile('TOVademecum','Interpolation',fileBase)
};

matFile = '';
for iFile = 1:numel(candidateFiles)
    if isfile(candidateFiles{iFile})
        matFile = candidateFiles{iFile};
        break
    end
end

if isempty(matFile)
    error(['Nao encontrei %s. Coloque este script na pasta do projeto ', ...
           'ou ajuste candidateFiles.'], fileBase);
end

fprintf('Arquivo carregado: %s\n', matFile);
S = load(matFile);

if ~isfield(S,'Interpolation')
    error('O arquivo nao contem a variavel Interpolation.');
end

I = S.Interpolation;

if ~isfield(I,'fun')
    error('Interpolation.fun nao foi encontrado.');
end

%% Faixas do teste
rhoMin = 1e-6;
rhoMax = 0.97;
nRho   = 301;
rhoVec = linspace(rhoMin,rhoMax,nRho).';

% Ajuste estes valores se o vademecum tiver outro intervalo em b.
bList = [-0.8, -0.4, 0, 0.4, 0.8];

% Grade mais densa para calcular a envolvente em b.
bEnvelope = linspace(min(bList),max(bList),161);

linearReference = rhoVec/rhoMax;

%% Componentes do tensor
components = struct( ...
    'name',  {'C1111','C2222','C1122','C1212'}, ...
    'index', {[1 1 1 1],[2 2 2 2],[1 1 2 2],[1 2 1 2]} );

results = struct();

%% Avaliacao e graficos por componente
for ic = 1:numel(components)

    compName = components(ic).name;
    idx      = components(ic).index;
    fun      = I.fun{idx(1),idx(2),idx(3),idx(4)};

    if isempty(fun)
        warning('%s esta vazio no vademecum. Componente ignorada.',compName);
        continue
    end

    Cabs  = zeros(nRho,numel(bList));
    Cnorm = zeros(nRho,numel(bList));
    ratio = nan(nRho,numel(bList));

    for ib = 1:numel(bList)
        bNow = bList(ib)*ones(nRho,1);
        Cnow = evaluateVademecumFunction(fun,bNow,rhoVec);

        Cabs(:,ib) = Cnow;

        Cref = evaluateVademecumFunction(fun,bList(ib),rhoMax);
        if abs(Cref) < 1e-14
            warning('%s: referencia quase nula para b = %.4g.',compName,bList(ib));
            Cnorm(:,ib) = NaN;
        else
            Cnorm(:,ib) = Cnow/Cref;
        end

        valid = linearReference > 1e-8;
        ratio(valid,ib) = Cnorm(valid,ib)./linearReference(valid);
    end

    results.(compName).Cabs  = Cabs;
    results.(compName).Cnorm = Cnorm;
    results.(compName).ratio = ratio;

    % Curvas normalizadas individualmente por C(b,rhoMax)
    figure('Name',[compName,' normalized']);
    hold on;
    for ib = 1:numel(bList)
        plot(rhoVec/rhoMax,Cnorm(:,ib), ...
            'LineWidth',1.4, ...
            'DisplayName',sprintf('b = %.1f',bList(ib)));
    end
    plot(rhoVec/rhoMax,linearReference,'k--','LineWidth',1.8, ...
        'DisplayName','linear: y = \rho/\rho_{max}');
    hold off;
    grid on;
    box on;
    xlabel('\rho/\rho_{max}');
    ylabel(sprintf('%s(b,\\rho)/%s(b,\\rho_{max})',compName,compName));
    title([compName,' - competitiveness of intermediate densities']);
    legend('Location','best');

    % Razao em relacao a resposta linear
    figure('Name',[compName,' competitiveness ratio']);
    hold on;
    for ib = 1:numel(bList)
        plot(rhoVec/rhoMax,ratio(:,ib), ...
            'LineWidth',1.4, ...
            'DisplayName',sprintf('b = %.1f',bList(ib)));
    end
    yline(1,'k--','LineWidth',1.8,'HandleVisibility','off');
    hold off;
    grid on;
    box on;
    xlabel('\rho/\rho_{max}');
    ylabel('normalized stiffness / linear reference');
    title([compName,' - ratio above/below linear response']);
    legend('Location','best');

    % Envolvente: melhor valor de b para cada rho.
    CenvAbs = zeros(nRho,numel(bEnvelope));
    for ib = 1:numel(bEnvelope)
        CenvAbs(:,ib) = evaluateVademecumFunction( ...
            fun,bEnvelope(ib)*ones(nRho,1),rhoVec);
    end

    % Normalizacao comum pela maior rigidez em rhoMax.
    CenvRef = max(CenvAbs(end,:));
    Cenv    = max(CenvAbs,[],2)/CenvRef;

    results.(compName).envelope = Cenv;

    figure('Name',[compName,' envelope']);
    plot(rhoVec/rhoMax,Cenv,'LineWidth',1.8, ...
        'DisplayName','max_b C(b,\rho)');
    hold on;
    plot(rhoVec/rhoMax,linearReference,'k--','LineWidth',1.8, ...
        'DisplayName','linear: y = \rho/\rho_{max}');
    hold off;
    grid on;
    box on;
    xlabel('\rho/\rho_{max}');
    ylabel(sprintf('max_b %s(b,\\rho) / max_b %s(b,\\rho_{max})', ...
        compName,compName));
    title([compName,' - envelope when b is free']);
    legend('Location','best');

end

%% Comparacao conjunta no caso b = 0
figure('Name','All components at b = 0');
hold on;

for ic = 1:numel(components)
    compName = components(ic).name;

    if ~isfield(results,compName)
        continue
    end

    [~,ib0] = min(abs(bList));
    plot(rhoVec/rhoMax,results.(compName).Cnorm(:,ib0), ...
        'LineWidth',1.5,'DisplayName',compName);
end

plot(rhoVec/rhoMax,linearReference,'k--','LineWidth',1.8, ...
    'DisplayName','linear: y = \rho/\rho_{max}');

hold off;
grid on;
box on;
xlabel('\rho/\rho_{max}');
ylabel('C_{ijkl}(0,\rho)/C_{ijkl}(0,\rho_{max})');
title('Normalized tensor components for b = 0');
legend('Location','best');

%% Resumo quantitativo
% Avalia as curvas em rho/rhoMax = 0.25, 0.50 e 0.75.
sampleFractions = [0.25 0.50 0.75];

fprintf('\n');
fprintf('============================================================\n');
fprintf('COMPETITIVENESS SUMMARY\n');
fprintf('ratio = [C/C(rhoMax)] / [rho/rhoMax]\n');
fprintf('ratio < 1: intermediate rho is below linear response\n');
fprintf('ratio = 1: approximately linear\n');
fprintf('ratio > 1: intermediate rho is above linear response\n');
fprintf('============================================================\n');

summaryRows = {};

for ic = 1:numel(components)
    compName = components(ic).name;

    if ~isfield(results,compName)
        continue
    end

    for ib = 1:numel(bList)
        for is = 1:numel(sampleFractions)
            [~,ir] = min(abs(rhoVec/rhoMax-sampleFractions(is)));

            summaryRows(end+1,:) = { ...
                compName, ...
                bList(ib), ...
                rhoVec(ir), ...
                results.(compName).Cnorm(ir,ib), ...
                results.(compName).ratio(ir,ib)}; %#ok<SAGROW>
        end
    end
end

summaryTable = cell2table(summaryRows, ...
    'VariableNames',{'Component','b','rho','NormalizedC','RatioToLinear'});

disp(summaryTable);

%% Testes de monotonicidade e curvatura
fprintf('\nDIAGNOSTICS\n');

for ic = 1:numel(components)
    compName = components(ic).name;

    if ~isfield(results,compName)
        continue
    end

    for ib = 1:numel(bList)
        y = results.(compName).Cnorm(:,ib);

        firstDiff  = diff(y);
        secondDiff = diff(y,2);

        nDecrease = nnz(firstDiff < -1e-8);
        fracBelow = mean(y(2:end-1) < linearReference(2:end-1)-1e-3);
        fracAbove = mean(y(2:end-1) > linearReference(2:end-1)+1e-3);

        fprintf(['%-6s, b=%5.2f: decreasing steps=%4d, ', ...
                 'below linear=%6.2f%%, above linear=%6.2f%%, ', ...
                 'mean second difference=% .3e\n'], ...
                 compName,bList(ib),nDecrease, ...
                 100*fracBelow,100*fracAbove,mean(secondDiff,'omitnan'));
    end
end

%% Salvar resultados
save('RhoCompetitiveness_Homogenizationtwovariables15.mat', ...
    'results','summaryTable','rhoVec','rhoMin','rhoMax', ...
    'bList','bEnvelope','components');

fprintf('\nResultados salvos em:\n');
fprintf('RhoCompetitiveness_Homogenizationtwovariables15.mat\n');

%% Funcao local robusta
function value = evaluateVademecumFunction(fun,b,rho)

    raw = fun(b,rho);

    % Interpolation.fun normalmente retorna um array numerico.
    % Este tratamento adicional evita erro caso alguma implementacao
    % retorne o valor dentro de uma cell.
    while iscell(raw)
        if isempty(raw)
            error('A funcao do vademecum retornou uma cell vazia.');
        end
        raw = raw{1};
    end

    value = full(raw);
    value = value(:);

    nExpected = max(numel(b),numel(rho));

    if numel(value) == 1 && nExpected > 1
        value = repmat(value,nExpected,1);
    end

    if numel(value) ~= nExpected
        error(['Tamanho inesperado retornado pelo vademecum: ', ...
               '%d valores para %d pontos.'],numel(value),nExpected);
    end
end
