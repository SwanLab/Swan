%% SCRIPT PER FER BATCHES DE SIMULACIONS

% TopOptTestTutorialDensityVolumePNorm
% TopOptTestTutorialDensityPerimeterPNorm
% TopOptTestTutorialLevelSetVolumePNorm
% TopOptTestTutorialLevelSetPerimeterPNorm

alpha   = [0.5 0.65 0.85 1];
pTarget = [0.5 0.65 0.85 1];

pNorm   = [1 4 8 16];

for ii = 1:size(pNorm,2)
    for jj = 1:size(alpha,2)
        p = pNorm(ii);
        a = alpha(jj);
        TopOptTestTutorialLevelSetVolumePNorm(p,a);
    end
end