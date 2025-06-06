%% SCRIPT PER FER BATCHES DE SIMULACIONS 1

values = [
        16 4 0.5 1
        ];
 
for ii = 1:size(values,1)
    p  = values(ii,1);
    pT = values(ii,2);
    gJ = values(ii,3);
    w  = values(ii,4);
    GripperProblemLevelSetPerimeterPNorm(p,pT,gJ,w);
end