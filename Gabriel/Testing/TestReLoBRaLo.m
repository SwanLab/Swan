clc
clear
close all

rng(1)

s.alpha = ...
    0.999;

s.temperature = ...
    0.1;

s.rhoExpected = ...
    0.5;

balancer = ...
    ReLoBRaLo(s);

lossValues = [...
    10.0    1.00;
    8.0    0.95;
    6.0    0.90;
    5.5    0.70;
    5.0    0.55;
    4.8    0.40];

nSteps = ...
    size(lossValues,1);

fprintf('\n')
fprintf('============================================\n')
fprintf('           RELOBRALO UNIT TEST\n')
fprintf('============================================\n')

for i = 1:nSteps

    losses = ...
        lossValues(i,:)';

    weights = ...
        balancer.update(losses);

    fprintf('\nStep %d\n',i)

    fprintf('Losses\n')
    disp(losses')

    fprintf('Weights\n')
    disp(weights')

    fprintf('Sum weights = %.15f\n',...
        sum(weights))

    if any(~isfinite(weights))

        error(...
            'Non-finite ReLoBRaLo weights.')

    end

    if any(weights < 0)

        error(...
            'Negative ReLoBRaLo weight.')

    end

    if abs(sum(weights)-length(weights)) > 1e-12

        error(...
            'ReLoBRaLo weights do not sum to nTerms.')

    end

end

fprintf('\n')
fprintf('Weight history\n')
disp(balancer.weightHistory)

fprintf('\n')
fprintf('Loss history\n')
disp(balancer.lossHistory)

fprintf('\n')
fprintf('Rho history\n')
disp(balancer.rhoHistory)

fprintf('\n')
fprintf('============================================\n')
fprintf('      RELOBRALO UNIT TEST PASSED\n')
fprintf('============================================\n')