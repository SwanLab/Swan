%% SCRIPT PER FER BATCHES DE SIMULACIONS

% TopOptTestTutorialDensityVolumePNorm
% TopOptTestTutorialDensityPerimeterPNorm
% TopOptTestTutorialLevelSetVolumePNorm
% TopOptTestTutorialLevelSetPerimeterPNorm

alpha   = [0.3 0.4 0.55];
pTarget = [0.3 0.4 0.5];

pNorm   = 16;

for ii = 1:size(pNorm,2)
    for jj = 1:size(pTarget,2)
        p  = pNorm(ii);
        pT = pTarget(jj);
        TopOptTestTutorialDensityPerimeterPNorm(p,pT);
    end
end