%% A

clc; clear; close all;

% r = 0.1;
r = 0.01:0.001:0.99;

aux2 = [1; 2; 3; 4; 5; 6; 7; 8];
S = [];


for j = 1:size(r,2)
    U = [];
    L = [];
    Object = {};
    auxl = [];
    
        % [~, u, l, mesh] = LevelSetInclusionAuto_raulHole(r(j),1);
        % string = strrep("UL_r"+r(j), ".", "_")+"-20x20-Hole"+".mat";

        [~, u, l, mesh] = LevelSetInclusionAuto_raul(r(j),1);
        string = strrep("UL_r"+num2str(r(j), '%.4f'), ".", "_")+"-20x20"+".mat";
        % 
        % [~, u, l, mesh] = LevelSetInclusionAuto_Ausetic();
        % string = "UL-Ausetic.mat";

        for i = 1:8
            auxl = cat(2, auxl, l(i, i:end));
        end
        disp(r(j))
        s = cat(2, r(j), auxl);
        S = cat(1, S, s);

        U         = cat(2, U, u);
        L         = cat(2, L, l);
        R         = r(j);
    
    
    
    save(string, "U", "L", "mesh","R");
    %save(string, "U", "L", "mesh");

end
writematrix(S, 'k.csv')