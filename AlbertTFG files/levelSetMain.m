%% A

clc; clear; close all; format long;


nPoints = 20;
distanceBetweenNodes = 2/nPoints;


%r = 0.1;
%r = 0.01:0.01:0.99;
r = 0.1:0.05:0.95;
%r = distanceBetweenNodes/2:distanceBetweenNodes/20:0.999;

sizeR = size(r, 2);

t1 = 'radius = ';
t2 = '\t';
t3 = 'Progress = ';
t4  = 'Time = ';
t5 = 's';

aux = ones(8,1);
identity = eye(8);
data2 = [];
data3 = [];
dataU = [];
dataPattern = [];
dataNormal = [];
time = 0;
timeInc = [];
tic
for j = 1:size(r,2)
    

        %% To save the data with the inclusion
        % [~, u, l, mesh] = LevelSetInclusionAuto_raulHole(r(j),1);
        % string = strrep("UL_r"+r(j), ".", "_")+"-20x20-Hole"+".mat";
        
        %% To save the data with the virtual inclusion
        
        [~, u, l, mesh] = LevelSetInclusionAuto_raul(r(j),nPoints);

        % dataU = cat(2, dataU, u);
        % string = strrep("UL_r"+num2str(r(j), '%.4f'), ".", "_")+"-20x20"+".mat";
        % string = strrep("r"+num2str(r(j), '%.2f'), ".", "_")+"-"+num2str(nPoints)+"x"+num2str(nPoints)+".mat";
        % U = u;
        % L = l;
        % R = r(j);
        % save(string, "U", "L", "mesh","R");
        
        % %save(string, "U", "L", "mesh");
        
        %% K network case 1
        % data = cat(2, identity, l);
        % string = "K_r"+num2str(r(j), '%.4f')+"-"+num2str(nPoints)+"x"+num2str(nPoints)+".csv";
        % writematrix(data, string)

        %% K network case 2
        inBetween = cat(2, r(j).*aux, identity, l);
        data2 = cat(1, data2, inBetween);

        %% K network case 3
        inBetween = r(j);
        % for i = 1:8
        %     inBetween = cat(2, inBetween, l(i,:));
        % end
        % data3 = cat(1, data3, inBetween);


        %% K network pattern recognition
        % patternRecognition = zeros(8,8);
        % 
        % currentNum = 1;
        % currentNumVal = max(max(l));
        % 
        % inBetween = cat(2, r(j).*aux, l);
        % data = cat(1, data, inBetween);
        % 
        % while max(max(ismember(patternRecognition, 0))) == 1
        %     for m = 1:8
        %         for n = 1:8
        %             if abs(abs(currentNumVal)-abs(l(m,n)))< 5*10^-3
        %                 patternRecognition(m,n) = currentNumVal/l(m,n)*currentNum;
        %             end
        %         end
        %     end
        % 
        %     looking = true;
        % 
        %     for m = 1:8
        %         for n = 1:8
        %             if patternRecognition(m,n)==0 && looking == true
        %                 currentNum = currentNum +1;
        %                 currentNumVal = l(m,n);
        %                 looking = false;
        %             end
        %         end
        %     end
        % 
        % 
        % end
        % 
        % inBetweenPattern = cat(2, r(j).*aux, patternRecognition);
        % dataPattern = cat(1, dataPattern, inBetweenPattern);
        % 
        % normalized = l./(max(abs(l)));
        % inBetweenNormal  = cat(2, r(j).*aux, normalized);
        % dataNormal = cat(1, dataNormal, inBetweenNormal);
    % 
    

     %% K network case 4
        % inBetween = r(j);
        %     inBetween = cat(2, inBetween, l(1,[1, 2, 3, 4, 5, 7]));
        % data = cat(1, data, inBetween);
    string = sprintf([t1, num2str(r(j)),t2,t3,num2str(j/sizeR),t2, t2, t4, num2str(toc), t5]);
    time = cat(2, time, toc);
    timeInc = cat(2, timeInc, toc-time(j) );
    disp(string);

end

%% K network case 2
string2 = "K_data case 2 20x20 v2 test.csv";

%% K network case 3
% string3 = "K_data case 3 20x20 v2 test.csv";

%% K network case 4
% string = "K_data case 4.csv";


%% Case 2; 3 and 4
writematrix(data2, string2)
% writematrix(data3, string3)

%% Pattern recognition
% writematrix(dataPattern, "PatternRecognition.csv")
% writematrix(dataNormal, "PatternNormalized.csv")
