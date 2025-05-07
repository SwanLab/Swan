%% A

clc; clear; close all;

% r = 0.1;
r = 0.05:0.05:0.95;


for j = 1:size(r,2)
    U = [];
    L = [];
    Object = {};
    
    %for i = 1:8
        % [~, u, l, mesh] = LevelSetInclusionAuto_raulHole(r(j),1);
        % string = strrep("UL_r"+r(j), ".", "_")+"-20x20-Hole"+".mat";

        % [~, u, l, mesh] = LevelSetInclusionAuto_raul(r(j),1);
        % string = strrep("UL_r"+r(j), ".", "_")+"-20x20"+".mat";
        % 
        [~, u, l, mesh] = LevelSetInclusionAuto_Ausetic();
        string = "UL-Ausetic.mat";

        U         = cat(2, U, u);
        L         = cat(2, L, l);
        R         = r(j);
    %end
    
    
    % save(string, "U", "L", "mesh","R");
    save(string, "U", "L", "mesh");

end