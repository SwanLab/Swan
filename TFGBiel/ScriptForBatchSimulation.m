%% SCRIPT PER FER BATCHES DE SIMULACIONS

% TopOptTestTutorialDensityVolumePNorm
% TopOptTestTutorialDensityPerimeterPNorm
% TopOptTestTutorialLevelSetVolumePNorm
% TopOptTestTutorialLevelSetPerimeterPNorm

alpha   = [0.5 0.65 0.85 1];
pTarget = [0.5 0.65 0.85];

pNorm   = [4 8];

for ii = 1:size(pNorm,2)
    for jj = 1:size(pTarget,2)
        p  = pNorm(ii);
        pt = pTarget(jj);
        TopOptTestTutorialDensityPerimeterPNorm(p,pt);
    end
end