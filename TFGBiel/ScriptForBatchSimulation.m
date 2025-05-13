%% SCRIPT PER FER BATCHES DE SIMULACIONS

% TopOptTestTutorialDensityVolumePNorm
% TopOptTestTutorialDensityPerimeterPNorm
% TopOptTestTutorialDensityIsoPerimetricPNorm
% TopOptTestTutorialLevelSetVolumePNorm
% TopOptTestTutorialLevelSetPerimeterPNorm
% TopOptTestTutorialLevelSetIsoPerimetricPNorm

% alpha   = [0.3 0.4 0.55];
% pTarget = [0.5 1 1.5 2 2.5];
% C       = [1.5 2 2.5 3 3.5 4 4.5];
% 
% pNorm   = [1 4 8 16];
% 
% for ii = 1:size(pNorm,2)
%     for jj = 1:size(C,2)
%         p  = pNorm(ii);
%         pT = C(jj);
%         TopOptTestTutorialDensityIsoPerimetricPNorm(p,C);
%     end
% end


values = [
            1  1
            4  2.5
            8  3
            16 4
            ];

for ii = 1:size(values,1)
    p = values(ii,1);
    C = values(ii,2);
    TopOptTestTutorialDensityIsoPerimetricPNorm(p,C);
end