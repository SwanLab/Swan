%% SCRIPT PER FER BATCHES DE SIMULACIONS

% TopOptTestTutorialDensityVolumePNorm
% TopOptTestTutorialDensityPerimeterPNorm
% TopOptTestTutorialLevelSetVolumePNorm
% TopOptTestTutorialLevelSetPerimeterPNorm

% alpha   = [0.3 0.4 0.55];
% pTarget = [0.5 1 1.5 2 2.5];
% 
% pNorm   = [1 4 8 16];

% for ii = 1:size(pNorm,2)
%     for jj = 1:size(pTarget,2)
%         p  = pNorm(ii);
%         pT = pTarget(jj);
%         TopOptTestTutorialLevelSetPerimeterPNorm(p,pT);
%     end
% end


values = [
            1 0.83
            4 1.9
            4 1.93
            8 2.3
            8 2.34
            8 2.4
            16 2.7
            16 2.74
            16 2.8
            ];

for ii = 1:size(values,1)
    p  = values(ii,1);
    pT = values(ii,2);
    TopOptTestTutorialLevelSetPerimeterPNorm(p,pT);
end