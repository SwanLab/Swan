%% A

clc; clear; close all;

r = [0.1, 0.3, 0.4];


for j = 1:size(r,2)
    U = [];
    L = [];
    
    %for i = 1:8
        [~, u, l] = LevelSetInclusionAuto_raulHole(r(j),1);
        string = strrep("UL_r"+r(j), ".", "_")+"-20x20-Hole"+".mat";

        % [~, u, l] = LevelSetInclusionAuto_raul(r(j),1);
        % string = strrep("UL_r"+r(j), ".", "_")+"-20x20"+".mat";

        %[~, u, l] = LevelSetInclusionAuto_Ausetic();
        %string = "UL-Ausetic.mat";

        U         = cat(2, U, u);
        L         = cat(2, L, l);
    %end
    
    
    save(string, "U", "L");
end