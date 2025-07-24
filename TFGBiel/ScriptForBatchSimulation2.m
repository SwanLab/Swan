%% SCRIPT PER FER BATCHES DE SIMULACIONS 2

values2 = [
        16 4.5 0.5 1
        ];
 
for ii = 1:size(values2,1)
    p  = values2(ii,1);
    C = values2(ii,2);
    gJ = values2(ii,3);
    w  = values2(ii,4);
    GripperProblemLevelSetIsoPerimetricPNorm(p,C,gJ,w);
end