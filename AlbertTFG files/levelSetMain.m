%% A

clc; clear; close all;

% r = 0.1;
r = 0.85;

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

        %Designa un nom per cada linea corresponent a un radi
        string = strrep("UL_r"+num2str(r(j), '%.4f'), ".", "_")+"-20x20"+".mat"; 

        
        % [~, u, l, mesh] = LevelSetInclusionAuto_Ausetic();
        % string = "UL-Ausetic.mat";

        for i = 1:8
            auxl = cat(2, auxl, l(i, i:end));   % concantena en horitz els 36 elem de la Kcoarse
        end
        disp(r(j))
        s = cat(2, r(j), auxl);   % concatena el radi+36 elem en horitzontal
        S = cat(1, S, s);         % concatena les files verticalment de cada s corresponent a un diferent radi

        U         = cat(2, U, u); % concatena els desplaçaments
        L         = cat(2, L, l); % concatena les L
        R         = r(j);
    
    
    
    save(string, "U", "L", "mesh","R");  %guarda el workspace per cert radi
    %save(string, "U", "L", "mesh");

end
writematrix(S, 'k.csv')  %guarda l'arxiu
