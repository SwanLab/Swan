%% SCRIPT PER FER BATCHES DE SIMULACIONS

% TopOptTestTutorialDensityVolumePNorm
% TopOptTestTutorialDensityPerimeterPNorm
% TopOptTestTutorialDensityIsoPerimetricPNorm
% TopOptTestTutorialLevelSetVolumePNorm
% TopOptTestTutorialLevelSetPerimeterPNorm
% TopOptTestTutorialLevelSetIsoPerimetricPNorm


% alpha   = [0.3 0.4 0.55];
% pTarget = [1.5 2 2.5 3 3.5];
% C       = [1.5 2 2.5 3 3.5 4 4.5];
% 
% pNorm   = 16;
% 
% for ii = 1:size(pNorm,2)
%     for jj = 1:size(C,2)
%         p  = pNorm(ii);
%         pT = pTarget(jj);
%         TopOptTestTutorialLevelSetPerimeterPNorm(p,pT);
%     end
% end


values = [
        16 2   0.05
        16 2.5 0.05
        16 3   0.05
        16 2   2
        16 2.5 2
        16 3   2
        ];

for ii = 1:size(values,1)
    p  = values(ii,1);
    pT = values(ii,2);
    gJ = values(ii,3);
    TopOptTestTutorialDensityPerimeterPNorm(p,pT,gJ);
end