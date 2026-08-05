%% TestCompetitivenessRho_b06_CommonReference.m
% Competitividade de rho para b = 0.6 usando referencia comum:
% C_ijkl(b = 0, rho = rhoMax).
%
% Compara:
%   C_ijkl(0.6,rho) / C_ijkl(0,rhoMax)
% com:
%   y = rho/rhoMax
%
% Esta normalizacao preserva o ganho/perda absoluta causado por b = 0.6.

clearvars;
close all;
clc;

%% Arquivo
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
    error('Nao encontrei %s.',fileBase);
end

fprintf('Arquivo carregado: %s\n',matFile);

S = load(matFile);

if ~isfield(S,'Interpolation') || ~isfield(S.Interpolation,'fun')
    error('Interpolation.fun nao foi encontrado.');
end

I = S.Interpolation;

%% Parametros
bFixed     = 0.6;
bReference = 0.0;

rhoMin = 1e-6;
rhoMax = 0.97;
nRho   = 301;

rhoVec = linspace(rhoMin,rhoMax,nRho).';
rhoN   = rhoVec/rhoMax;

components = struct( ...
    'name',  {'C1111','C2222','C1122','C1212'}, ...
    'index', {[1 1 1 1],[2 2 2 2],[1 1 2 2],[1 2 1 2]} );

results = struct();

%% Figura conjunta com referencia comum
figure('Name','Normalized tensor components for b = 0.6 - common reference');
hold on;

for ic = 1:numel(components)

    name = components(ic).name;
    idx  = components(ic).index;

    fun = I.fun{idx(1),idx(2),idx(3),idx(4)};

    if isempty(fun)
        warning('%s esta vazio no vademecum.',name);
        continue
    end

    bVec = bFixed*ones(size(rhoVec));

    C = evaluateVademecumFunction(fun,bVec,rhoVec);

    % Referencia comum: mesma componente em b = 0 e rho = rhoMax
    Cref = evaluateVademecumFunction(fun,bReference,rhoMax);

    if abs(Cref) < 1e-14
        warning('%s possui referencia quase nula.',name);
        Cnorm = nan(size(C));
    else
        Cnorm = C/Cref;
    end

    ratio = nan(size(Cnorm));
    valid = rhoN > 1e-8;
    ratio(valid) = Cnorm(valid)./rhoN(valid);

    results.(name).C     = C;
    results.(name).Cref  = Cref;
    results.(name).Cnorm = Cnorm;
    results.(name).ratio = ratio;

    plot(rhoN,Cnorm,'LineWidth',1.6,'DisplayName',name);
end

plot(rhoN,rhoN,'k--','LineWidth',1.8, ...
    'DisplayName','linear: y = \rho/\rho_{max}');

hold off;
grid on;
box on;

xlabel('\rho/\rho_{max}');
ylabel('C_{ijkl}(0.6,\rho)/C_{ijkl}(0,\rho_{max})');
title('Tensor components for b = 0.6 with common reference b = 0');
legend('Location','best');

%% Razao de competitividade com referencia comum
figure('Name','Competitiveness ratio for b = 0.6 - common reference');
hold on;

for ic = 1:numel(components)

    name = components(ic).name;

    if ~isfield(results,name)
        continue
    end

    plot(rhoN,results.(name).ratio, ...
        'LineWidth',1.6, ...
        'DisplayName',name);
end

yline(1,'k--','LineWidth',1.8, ...
    'DisplayName','linear reference');

hold off;
grid on;
box on;

xlabel('\rho/\rho_{max}');
ylabel('[C(0.6,\rho)/C(0,\rho_{max})]/[\rho/\rho_{max}]');
title('Competitiveness ratio for b = 0.6 with common reference');
legend('Location','best');

%% Graficos individuais
for ic = 1:numel(components)

    name = components(ic).name;

    if ~isfield(results,name)
        continue
    end

    figure('Name',[name,' - b = 0.6 - common reference']);

    plot(rhoN,results.(name).Cnorm, ...
        'LineWidth',1.8, ...
        'DisplayName',sprintf('%s, b = %.1f',name,bFixed));

    hold on;

    plot(rhoN,rhoN,'k--','LineWidth',1.8, ...
        'DisplayName','linear: y = \rho/\rho_{max}');

    hold off;
    grid on;
    box on;

    xlabel('\rho/\rho_{max}');
    ylabel(sprintf('%s(0.6,\\rho)/%s(0,\\rho_{max})',name,name));
    title(sprintf('%s - b = 0.6 with common reference',name));
    legend('Location','best');
end

%% Resumo quantitativo
sampleFractions = [0.25 0.50 0.75];
summaryRows = {};

for ic = 1:numel(components)

    name = components(ic).name;

    if ~isfield(results,name)
        continue
    end

    for is = 1:numel(sampleFractions)

        [~,ir] = min(abs(rhoN-sampleFractions(is)));

        summaryRows(end+1,:) = { ...
            name, ...
            bFixed, ...
            bReference, ...
            rhoN(ir), ...
            rhoVec(ir), ...
            results.(name).Cnorm(ir), ...
            results.(name).ratio(ir)}; %#ok<SAGROW>
    end
end

summaryTable = cell2table(summaryRows, ...
    'VariableNames',{ ...
        'Component', ...
        'b', ...
        'bReference', ...
        'rhoNormalized', ...
        'rho', ...
        'NormalizedC', ...
        'RatioToLinear'});

fprintf('\n============================================================\n');
fprintf('COMPETITIVENESS FOR b = %.2f\n',bFixed);
fprintf('COMMON REFERENCE: C_ijkl(b = %.2f, rho = %.2f)\n', ...
    bReference,rhoMax);
fprintf('Ratio < 1: below linear response\n');
fprintf('Ratio = 1: linear response\n');
fprintf('Ratio > 1: above linear response\n');
fprintf('============================================================\n');

disp(summaryTable);

%% Diagnostico
fprintf('\nDIAGNOSTICS FOR b = %.2f WITH COMMON REFERENCE\n',bFixed);

for ic = 1:numel(components)

    name = components(ic).name;

    if ~isfield(results,name)
        continue
    end

    y = results.(name).Cnorm;

    fracBelow = mean(y(2:end-1) < rhoN(2:end-1)-1e-3);
    fracAbove = mean(y(2:end-1) > rhoN(2:end-1)+1e-3);

    fprintf('%-6s: below linear = %6.2f%%, above linear = %6.2f%%\n', ...
        name,100*fracBelow,100*fracAbove);
end

%% Salvar
save('RhoCompetitiveness_b06_CommonReference.mat', ...
    'results', ...
    'summaryTable', ...
    'rhoVec', ...
    'rhoN', ...
    'rhoMin', ...
    'rhoMax', ...
    'bFixed', ...
    'bReference');

fprintf('\nResultados salvos em:\n');
fprintf('RhoCompetitiveness_b06_CommonReference.mat\n');

%% Funcao local
function value = evaluateVademecumFunction(fun,b,rho)

    raw = fun(b,rho);

    while iscell(raw)
        if isempty(raw)
            error('A funcao retornou uma cell vazia.');
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
               '%d valores para %d pontos.'], ...
               numel(value),nExpected);
    end
end
