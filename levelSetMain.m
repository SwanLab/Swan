%% A

clc; clear; close all;

r = [0.1, 0.3, 0.4];


for j = 1:size(r,2)
    U = [];
    L = [];
    
    for i = 1:8
        [~, u, l] = LevelSetInclusionAuto(r(j), i);
        U         = cat(2, U, u);
        L         = cat(2, L, l);
    end
    string = strrep("UL_r"+r(j), ".", "_")+".mat";
    save(string, "U", "L");
end